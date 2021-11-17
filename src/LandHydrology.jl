module LandHydrology
using CLIMAParameters
using DocStringExtensions
using ClimaCore:Fields
abstract type AbstractLandModel end

# Previously, all functions in SubComponentModels applied to any AbstractModel, and LandHydrologyModel <: AbstractModel. Aside from
# consistency, what is the benefit of this ? This would make a lot of sense if we could run each model independently of the LandModel
# wrapper. But I'm not sure we can. 
include("SubComponentModels.jl")
using .SubComponentModels
include("Domains/Domains.jl")
export make_rhs, LandHydrologyModel, PrescribedAtmosState, NoAtmosState, initialize_land_states

"""
    AbstractAtmosState

Abstract type encompassing land-atmosphere interactions.
"""
abstract type AbstractAtmosState end

"""
     PrescribedAtmosState{FT<:AbstractFloat, p} <: AbstractAtmosState

The concrete type to be chosen when the land model is driven by prescribed atmospheric variables,
which can (eventually) be functions of space and time. 
"""
struct PrescribedAtmosState{FT<:AbstractFloat, p} <: AbstractAtmosState
    precipitation::p
    PrescribedAtmosState{FT}(precipitation::p) where {FT, p} = new{FT,p}(precipitation)
end

"""
    CoupledAtmosState <: AbstractAtmosState

The concrete type to be chosen when coupling the land model to a dynamic
atmosphere model.
"""
struct CoupledAtmosState <: AbstractAtmosState end

"""
    NoAtmosState <: AbstractAtmosState

A concrete type to be chosen when no interactions with the atmosphere are required to simulate
the land surface model.
"""
struct NoAtmosState <: AbstractAtmosState end

"""
    LandHydrologyModel{FT, sm, sfcm, ep, atm} <: AbstractLandModel

The concrete type which defines all of the physics and setup of the land 
hydrology model, including domains, parameters, prognostic variables, differential equations, and 
included interactions between subcomponents.
# Fields
$(DocStringExtensions.FIELDS)

"""
struct LandHydrologyModel{FT<: AbstractFloat,
                          sm <: AbstractModel,
                          sfcm<:AbstractModel,
                          ep <:AbstractEarthParameterSet,
                          atm <:AbstractAtmosState
                          } <: AbstractLandModel
    "The soil model"
    soil::sm
    "The surface water model"
    sfc_water::sfcm
    "Earth Parameter Set"
    earth_param_set::ep
    "Atmospheric state - either coupled, prescribed, or unused"
    atmos_state::atm
    LandHydrologyModel{FT}(soil::sm,
                           sfc_water::sfcm,
                           earth_param_set::ep;
                           atmos_state::atm = NoAtmosState()
                           ) where {FT, sm, sfcm,ep,atm} = new{FT, sm, sfcm, ep, atm}(
                               soil, sfc_water,earth_param_set,atmos_state)

end

"""
    function make_land_update_aux(model::LandHydrologyModel)

Creates the update_aux!(Ya, Y, t) function for the entire land model, which is called
at the beginning of each rhs calculation.

Auxilary variables are functions of time, space, and the prognostic variables.
They can used to incorporate prescribed parameters (driving the system), to store
functions of state that are required often in RHS calculations, or to solve for
variables that are updated via algebraic equations. Each time the right hand side
function is computed, the auxiliary variables must be updated first so that their
values correspond to the current time `t`. 
"""
function make_land_update_aux(model::LandHydrologyModel)
    interactions_update_aux! = SubComponentModels.make_update_aux(model.atmos_state, model.soil, model.sfc_water, model)
    soil_update_aux! = SubComponentModels.make_update_aux(model.soil, model)
    sfc_water_update_aux! = SubComponentModels.make_update_aux(model.sfc_water, model)
    function update_aux!(Ya, Y, t)
        interactions_update_aux!(Ya, Y, t)
        soil_update_aux!(Ya, Y, t)
        sfc_water_update_aux!(Ya,Y, t)
    end
    return update_aux!
end


"""
    function make_rhs(model::LandHydrologyModel)

Creates the rhs!(dY, Y, Ya t) function for the entire land model,
where dY is the rhs of the ODE evaluated at the current prognostic state `Y`,
time `t`, and parameters `Ya`. 

That is, if your system of equations is

``
∂Y_i/∂t = dY_i = f(Y, Ya(x,y,z,t,Y), t),
``

rhs! is a function which updates `dY` in place with the time derivative of the
variables `Y` at the current time `t`.  This requires updating the parameter `Y` 
so that their values correspond to the current time, which is why the rhs! function
for the `LandHydrologyModel` invokes `update_aux!` prior to updating `dY`.

If your ODEs result from the semi-discrete form of a PDE,
the rhs function must give the rhs of the discretized system. 
"""
function make_rhs(model::LandHydrologyModel)
    update_aux! = make_land_update_aux(model)
    soil_tendency_terms! = SubComponentModels.make_tendency_terms(model.soil, model)
    sfc_water_tendency_terms! = SubComponentModels.make_tendency_terms(model.sfc_water, model)
    function rhs!(dY, Y, Ya, t)
        update_aux!(Ya,Y,t)
        soil_tendency_terms!(dY,Y, Ya, t)
        sfc_water_tendency_terms!(dY,Y,Ya,t)
    end
    return rhs!
end


"""
    initialize_land_states(model::LandHydrologyModel, f::NamedTuple)

Initialize the land model's prognostic and auxiliary states, by calling
the `SubComponentModels.initialize_states` function for each subcomponent.

The prognostic variables are initialized according to the functions passed
in via the NamedTupled `f`, which is expected to have keys that match
the names of the subcomponent models.

The auxiliary parameters are initalized with all zero values. In the computation
of the `rhs` on the first step, these will be updated with their values at the
initial time.
"""
function initialize_land_states(model::LandHydrologyModel, f::NamedTuple)
    subcomponents = (:soil, :sfc_water)
    Y = Dict()
    Ya = Dict()
    for sc_name in subcomponents
        sc_model = getproperty(model, sc_name)
        if typeof(sc_model) != NotIncluded
            Y_sc, Ya_sc = SubComponentModels.initialize_states(sc_model, model, getproperty(f, sc_name))
            if sizeof(Y_sc) > 0.0
                push!(Y, sc_name => Y_sc)
            end
            if sizeof(Ya_sc) >0.0
                push!(Ya, sc_name => Ya_sc)
            end
        end
        
        
    end
    # do we always want this? should we only add this for certain types?
    push!(Ya, :soil_infiltration => 0.0,)
    return Fields.FieldVector(; Y...),Fields.FieldVector(; Ya...)
end

"""
    default_land_initial_conditions(model::LandHydrologyModel)

Initialize the land model's prognostic and auxiliary states with default
values, by calling the `SubComponentModels.default_initial_conditions`
 function for each subcomponent.
"""
function default_land_initial_conditions(model::LandHydrologyModel)
    subcomponents = (:soil, :sfc_water)
    Y = Dict()
    Ya = Dict()
    for sc_name in subcomponents
        sc_model = getproperty(model, sc_name)
        if typeof(sc_model) != NotIncluded
            Y_sc, Ya_sc = SubComponentModels.default_initial_conditions(sc_model, model)
            if sizeof(Y_sc) > 0.0
                push!(Y, sc_name => Y_sc)
            end
            if sizeof(Ya_sc) >0.0
                push!(Ya, sc_name => Ya_sc)
            end
        end
        
        
    end
    # do we always want this? should we only add this for certain types?
    push!(Ya, :soil_infiltration => 0.0,)
    return Fields.FieldVector(; Y...),Fields.FieldVector(; Ya...)
end


include(joinpath("SoilModel", "SoilInterface.jl"))
include(joinpath("SurfaceFlowModel", "SurfaceWater.jl"))
include("Interactions/Interactions.jl")
include("Simulations/Simulations.jl")

end # module
