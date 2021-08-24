export SoilHydrologyModel,
    SoilModel, PrescribedTemperatureModel, SoilEnergyModel

abstract type AbstractSoilModel end

struct SoilEnergyModel <: AbstractSoilModel
    initial_conditions::Any
    boundary_conditions::Any
end
struct SoilHydrologyModel <: AbstractSoilModel
    initial_conditions::Any
    boundary_conditions::Any
end

struct PrescribedTemperatureModel <: AbstractSoilModel
    profile::Any
end


struct SoilModel{A, B} <: AbstractSoilModel
    energy_model::AbstractSoilModel
    hydrology_model::AbstractSoilModel
    soil_params::A
    earth_params::B
end
