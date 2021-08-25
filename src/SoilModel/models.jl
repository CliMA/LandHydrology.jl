export SoilHydrologyModel,
    SoilModel,
    PrescribedTemperatureModel,
    PrescribedHydrologyModel,
    SoilEnergyModel

abstract type AbstractSoilModel end
"""
    SoilEnergyModel{ic,bc}

The model type to be used when the user wants to simulate
heat transfer in soil by solving the heat partial differential equation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilEnergyModel{ic, bc} <: AbstractSoilModel
    "function returning initial conditions for `ρe_int`"
    initial_conditions::ic
    "Boundary conditions for `ρe_int`"
    boundary_conditions::bc
end

"""
    SoilHydrologyModel{ic, bc}

The model type to be used when the user wants to simulate
the flow of water in soil by solving Richards equation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilHydrologyModel{ic, bc} <: AbstractSoilModel
    "function returning initial conditions for `ϑ_l` and `θ_i`"
    initial_conditions::ic
    "Boundary conditions for `ϑ_l`"
    boundary_conditions::bc
end

"""
    Base.@kwdef struct PrescribedTemperatureModel <: AbstractSoilModel
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
Base.@kwdef struct PrescribedTemperatureModel <: AbstractSoilModel
    "Profile of (z,t) for temperature"
    T_profile::Function = (z,t) -> eltype(z)(288)
end



"""
    Base.@kwdef struct PrescribedHydrologyModel <: AbstractSoilModel
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
Base.@kwdef struct PrescribedHydrologyModel <: AbstractSoilModel
    "Profile of (z,t) for ϑ_l"
    ϑ_l_profile::Function = (z,t) -> eltype(z)(0.0)
    "Profile of (z,t) for θ_i"
    θ_i_profile::Function = (z,t) -> eltype(z)(0.0)
end



"""
    SoilModel{em <: AbstractSoilModel, hm <: AbstractSoilModel, A,B}

The model type for the soil model. This is a composition of the individual
types, i.e. the energy model type and hydrology model type, as well as the
parameters needed by the soil model.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilModel{em <: AbstractSoilModel, hm <: AbstractSoilModel, A, B} <:
       AbstractSoilModel
    "Soil energy model - prescribed or dynamics"
    energy_model::em
    "Soil hydrology model - prescribed or dynamic"
    hydrology_model::hm
    "Soil parameters"
    soil_params::A
    "Earth parameter set"
    earth_params::B
end
