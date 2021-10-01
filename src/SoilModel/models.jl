export SoilHydrologyModel,
    SoilModel,
    PrescribedTemperatureModel,
    PrescribedHydrologyModel,
    SoilEnergyModel

abstract type AbstractSoilComponentModel end
"""
    SoilEnergyModel
The model type to be used when the user wants to simulate
heat transfer in soil by solving the heat partial differential equation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilEnergyModel <: AbstractSoilComponentModel end

"""
    SoilHydrologyModel

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
    Base.@kwdef struct PrescribedTemperatureModel <: AbstractSoilComponentModel
        "Profile of (z,t) for temperature"
         T_profile::Function = (z,t) -> eltype(z)(288)
    end
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
    Base.@kwdef struct PrescribedHydrologyModel <: AbstractSoilComponentModel
        "Profile of (z,t) for ϑ_l"
        ϑ_l_profile::Function = (z,t) -> eltype(z)(0.0)
        "Profile of (z,t) for θ_i"
        θ_i_profile::::Function = (z,t) -> eltype(z)(0.0)
    end
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
    SoilModel{domain, em <: AbstractSoilModel, hm <: AbstractSoilModel, bc, A,B}

The model type for the soil model.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SoilModel{
    dm,
    em <: AbstractSoilComponentModel,
    hm <: AbstractSoilComponentModel,
    bc,
    A,
    B,
} <: AbstractModel
    domain::dm
    "Soil energy model - prescribed or dynamics"
    energy_model::em
    "Soil hydrology model - prescribed or dynamic"
    hydrology_model::hm
    "Boundary conditions tuple"
    boundary_conditions::bc
    "Tuple of sources eventually"
    sources::AbstractLandSource = NoSource()
    "Soil parameters"
    soil_param_set::A
    "Earth parameter set"
    earth_param_set::B
    "name"
    name::Symbol = :soil
    "variables"
    variables::Tuple = (:ϑ_l, :θ_i, :ρe_int)
end


