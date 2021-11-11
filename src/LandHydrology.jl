# Deal with default_initial_conditions
module LandHydrology
using CLIMAParameters
using DocStringExtensions
using ClimaCore:Fields
struct EarthParameterSet <: AbstractEarthParameterSet end

include("Domains/Domains.jl")
include("Models.jl")
using .Models: AbstractModel, NotIncluded
import .Models:  make_tendency_terms, default_initial_conditions, make_update_aux, initialize_states
export make_rhs, LandHydrologyModel, PrescribedAtmosState, NoAtmosState

abstract type AbstractAtmosState end

struct PrescribedAtmosState{FT<:AbstractFloat, p} <: AbstractAtmosState
    precipitation::p
    PrescribedAtmosState{FT}(precipitation::p) where {FT, p} = new{FT,p}(precipitation)
end

struct CoupledAtmosState <: AbstractAtmosState end # signal to be coupled to dynamic model
struct NoAtmosState <: AbstractAtmosState end # signal to be coupled to dynamic model

"""
    LandHydrologyModel{FT<: AbstractFloat, sm <: AbstractModel} <: AbstractModel

The structure which defines all of the physics and setup of the land 
hydrology model, including domains, parameters, prognostic variables, differential equations, and included interactions between subcomponents.
# Fields
$(DocStringExtensions.FIELDS)

"""
struct LandHydrologyModel{FT<: AbstractFloat, sm <: AbstractModel, sfcm<:AbstractModel, atm <:AbstractAtmosState} <: AbstractModel
    "The soil model"
    soil::sm
    "The surface water model"
    sfc_water::sfcm
    "Atmospheric state - either coupled, prescribed, or unused"
    atmos_state::atm
    LandHydrologyModel{FT}(soil::sm, sfc_water::sfcm; atmos_state::atm = NoAtmosState()) where {FT, sm, sfcm,atm} = new{FT, sm, sfcm, atm}(
    soil, sfc_water,atmos_state)
end

"""
    function Models.make_update_aux(model::LandHydrologyModel)

Creates the update_aux!(Ya, Y, t) function for the entire land model, which is called
at the beginning of each rhs calculation.

Auxilary variables are functions of time, space, and the prognostic variables.
They can used to incorporate prescribed parameters (driving the system), to store
functions of state that are required often in RHS calculations, or to solve for
variables that are updated via algebraic equations. Each time the right hand side
function is computed, the auxiliary variables must be updated first so that their
values correspond to the current time `t`. 
"""
function Models.make_update_aux(model::LandHydrologyModel)
    soil_update_aux! = Models.make_update_aux(model.soil, model)
    sfc_water_update_aux! = Models.make_update_aux(model.sfc_water, model)
    function update_aux!(Ya, Y, t)
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

If you do not define a 
method of `make_rhs!` for your specific model type, the default is not
to update the prognostic state `Y`: `dY = 0`.
"""
function make_rhs(model::LandHydrologyModel)
    update_aux! = make_update_aux(model)
    soil_tendency_terms! = Models.make_tendency_terms(model.soil, model)
    sfc_water_tendency_terms! = Models.make_tendency_terms(model.sfc_water, model)
    function rhs!(dY, Y, Ya, t)
        update_aux!(Ya,Y,t)
        soil_tendency_terms!(dY,Y, Ya, t)
        sfc_water_tendency_terms!(dY,Y,Ya,t)
    end
    return rhs!
end



function Models.initialize_states(model::LandHydrologyModel, f::NamedTuple, t0::Real)
    subcomponents = (:soil, :sfc_water)
    Y = Dict()
    Ya = Dict()
    for sc_name in subcomponents
        sc_model = getproperty(model, sc_name)
        if typeof(sc_model) != NotIncluded
            Y_sc, Ya_sc = Models.initialize_states(sc_model, getproperty(f, sc_name), t0)
            push!(Y, sc_name => getproperty(Y_sc, sc_name))
            push!(Ya, sc_name => getproperty(Ya_sc, sc_name))
        end
        
    end
    return Fields.FieldVector(; Y...),Fields.FieldVector(; Ya...)
end

function Models.default_initial_conditions(model::LandHydrologyModel)
    subcomponents = (:soil, :sfc_water)
    Y = Dict()
    Ya = Dict()
    for sc_name in subcomponents
        sc_model = getproperty(model, sc_name)
        sc_model = getproperty(model, sc_name)
        if typeof(sc_model) != NotIncluded
            Y_sc, Ya_sc = Models.default_initial_conditions(sc_model)
            push!(Y, sc_name => getproperty(Y_sc, sc_name))
            push!(Ya, sc_name => getproperty(Ya_sc, sc_name))
        end
        
    end
    return Fields.FieldVector(; Y...),Fields.FieldVector(; Ya...)
end



include(joinpath("SoilModel", "SoilInterface.jl"))
include(joinpath("SurfaceFlowModel", "SurfaceWater.jl"))
include("Simulations/Simulations.jl")

end # module
