export RichardsEquation,SoilHeatEquation, SoilModel, NoSoilModel

abstract type AbstractSoilModel end

"""
    struct SoilHeatEquation{bc, a, b, lp, ip} <: AbstractSoilModel
        "Boundary conditions tuple"
        boundary_conditions::bc
        "Soil parameters"
        soil_param_set::a
        "Earth parameter set"
        earth_param_set::b
        "Profile of (z,t) for θ_l"
        θ_l_profile::lp
        "Profile of (z,t) for θ_i"
        θ_i_profile::ip
    end

The model type to be used when the user does not wish to solve
the Richards equation, but to solve the PDE for heat transfer in soil.

In this case, the user may specify a water content profile in the soil. 
The default is completely dry soil.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilHeatEquation{bc, a, b, lp, ip} <: AbstractSoilModel
    "Boundary conditions tuple"
    boundary_conditions::bc
    "Soil parameters"
    soil_param_set::a
    "Earth parameter set"
    earth_param_set::b
    "Profile of (z,t) for θ_l"
    θ_l_profile::lp
    "Profile of (z,t) for θ_i"
    θ_i_profile::ip
end

function SoilHeatEquation(
    boundary_conditions,
    soil_param_set,
    earth_param_set;
    θ_l_profile::Function = (z,t) -> eltype(z)(0.0),
    θ_i_profile::Function = (z,t) -> eltype(z)(0.0),
)
    args = (boundary_conditions, soil_param_set, earth_param_set, θ_l_profile, θ_i_profile)
    return SoilHeatEquation{typeof.(args)...}(args...)
end


"""
    struct RichardsEquation{bc, a, b, tp} <: AbstractSoilModel
        "Boundary conditions tuple"
        boundary_conditions::bc
        "Soil parameters"
        soil_param_set::a
        "Earth parameter set"
        earth_param_set::b
        "Profile of (z,t) for T"
        T_profile::tp
    end

The model type to be used when the user wants to solve Richards equation.

The user may optionally supply a temperature profile in the soil in order to
mode the temperature dependence of conductivity. The default is no effect.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct RichardsEquation{bc, a, b, tp} <: AbstractSoilModel
    "Boundary conditions tuple"
    boundary_conditions::bc
    "Soil parameters"
    soil_param_set::a
    "Earth parameter set"
    earth_param_set::b
    "Profile of (z,t) for T"
    T_profile::tp
end


function RichardsEquation(
    boundary_conditions,
    soil_param_set,
    earth_param_set;
    T_profile::Function = (z,t) -> eltype(z)(288.0)
)
    args = (boundary_conditions, soil_param_set, earth_param_set, T_profile)
    return RichardsEquation{typeof.(args)...}(args...)
end


"""
    NoSoilModel <: AbstractSoilModel

To be used when no dynamic soil model is required.
"""
struct NoSoilModel <: AbstractSoilModel end

"""
    SoilModel{bc, A,B}

The model type for the soil model, when both energy and water flow are simulated prognostically.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilModel{bc, A, B} <:
       AbstractSoilModel
    "Boundary conditions tuple"
    boundary_conditions::bc
    "Soil parameters"
    soil_param_set::A
    "Earth parameter set"
    earth_param_set::B
end


include("boundary_conditions.jl")
include("right_hand_side.jl")
