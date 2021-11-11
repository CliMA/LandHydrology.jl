export SoilHydrologyModel,
    SoilModel,
    PrescribedTemperatureModel,
    PrescribedHydrologyModel,
    SoilEnergyModel,
    compute_infiltration

function compute_infiltration(_,_,_) return 0.0 end


abstract type AbstractSoilComponentModel end
"""
    SoilEnergyModel <: AbstractSoilComponentModel

The model type to be used when the user wants to simulate
heat transfer in soil by solving the heat partial differential equation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilEnergyModel <: AbstractSoilComponentModel end

"""
    SoilHydrologyModel <: AbstractSoilComponentModel

The model type to be used when the user wants to simulate
the flow of water in soil by solving Richards equation.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SoilHydrologyModel{FT <: AbstractFloat} <:
                   AbstractSoilComponentModel
    hydraulic_model::vanGenuchten{FT} = vanGenuchten{FT}()
    viscosity_factor::AbstractConductivityFactor{FT} = NoEffect{FT}()
    impedance_factor::AbstractConductivityFactor{FT} = NoEffect{FT}()
end


"""
    PrescribedTemperatureModel <: AbstractSoilComponentModel

The model type to be used when the user does not wish to solve
the heat partial differential equation, but instead wishes to prescibe
a temperature profile in the soil.

This is useful for situations where Richards Equation alone is sufficient.
Because the hydraulic conductivity can be a function of temperature, a
temperature profile can be supplied in order to simulate that.
The default is 288K across the domain, the reference temperature for the viscosity effect.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PrescribedTemperatureModel <: AbstractSoilComponentModel
    "Profile of (z,t) for temperature"
    T_profile::Function = (z, t) -> eltype(t)(288)
end



"""
    PrescribedHydrologyModel <: AbstractSoilComponentModel

The model type to be used when the user does not wish to solve
the Richards equation, but instead wishes to prescibe
a water content profile in the soil.

This is useful for situations where only the heat equation is to be solved.
Because the thermal conductivity and heat capacities depend on water content,
a water profile must be defined to solve the heat equation.
The default for both ice
and liquid water content is zero, which applies for totally dry soil. 
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PrescribedHydrologyModel <: AbstractSoilComponentModel
    "Profile of (z,t) for ϑ_l"
    ϑ_l_profile::Function = (z, t) -> eltype(z)(0.0)
    "Profile of (z,t) for θ_i"
    θ_i_profile::Function = (z, t) -> eltype(z)(0.0)
end



"""
    SoilModel{FT, domain, em <: AbstractSoilModel, hm <: AbstractSoilModel, bc, sp,ep,n}

The model type for the soil model.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilModel{FT, dm, em, hm, bc, sp, ep, n} <: AbstractModel
    domain::dm
    "Soil energy model - prescribed or dynamics"
    energy_model::em
    "Soil hydrology model - prescribed or dynamic"
    hydrology_model::hm
    "Boundary conditions tuple"
    boundary_conditions::bc
    "Soil parameters"
    soil_param_set::sp
    "Earth parameter set"
    earth_param_set::ep
    "name"
    name::n
end

function SoilModel(
    ::Type{FT};
    domain::AbstractVerticalDomain{FT},
    energy_model::AbstractSoilComponentModel,
    hydrology_model::AbstractSoilComponentModel,
    boundary_conditions::BC,
    soil_param_set::SoilParams{FT} = SoilParams{FT}(),
    earth_param_set::EarthParameterSet = EarthParameterSet(),
    name::Symbol = :soil,
) where {FT, BC}
    args = (
        domain,
        energy_model,
        hydrology_model,
        boundary_conditions,
        soil_param_set,
        earth_param_set,
        name,
    )
    return SoilModel{FT, typeof.(args)...}(args...)
end



"""
    default_initial_conditions(model::SoilModel{f, dm, SoilEnergyModel{f}, SoilHydrologyModel{f}})

Returns a default set of initial conditions for the soil.

The default is an isothermal soil, at the reference temperature T0,
no ice, and a volumetric water fraction constant throughout the domain,
at half of porosity.
"""
function Models.default_initial_conditions(
    model::SoilModel{f, dm, SoilEnergyModel, SoilHydrologyModel{f}},
) where {f, dm}
    function ic(z::f, m::SoilModel)
        param_set = model.earth_param_set
        T = f(273.16)
        θ_i = f(0.0)
        θ_l = f(0.5) * model.soil_param_set.ν
        ρcds = model.soil_param_set.ρc_ds
        ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρcds, param_set)
        ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, param_set)
        return (ϑ_l = θ_l, θ_i = θ_i, ρe_int = ρe_int)
    end
    t0 = f(0.0)
    return Models.initialize_states(model, ic, t0)
end

function Models.default_initial_conditions(model::SoilModel)
    error("No default IC exist for this type of soil model.")
end




function get_temperature(model::SoilModel{f, dm, SoilEnergyModel,}, Y::Fields.FieldVector, Ya::Fields.FieldVector) where {dm, f}
    ρe_int = Y.soil.ρe_int
    ϑ_l, θ_i = get_water_content(model.hydrology_model, Y, Ya)
    # Parameters
    sp = model.soil_param_set
    param_set = model.earth_param_set
    @unpack ν, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp
    
    # Compute center values of everything
    ν_eff = ν .- θ_i
    θ_l = volumetric_liquid_fraction.(ϑ_l, ν_eff)
    
    ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
    T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))
    return T
end


function get_temperature(model::SoilModel{f, dm, PrescribedTemperatureModel,}, Y::Fields.FieldVector, Ya::Fields.FieldVector) where {dm, f}
    return Ya.soil.T
end


function get_water_content(model::SoilHydrologyModel, Y::Fields.FieldVector,Ya::Fields.FieldVector)
    return Y.soil.ϑ_l, Y.soil.θ_i
end


function get_water_content(model::PrescribedHydrologyModel, Y::Fields.FieldVector,Ya::Fields.FieldVector)
    return Ya.soil.ϑ_l, Ya.soil.θ_i
end

include("boundary_conditions.jl")
include("right_hand_side.jl")
