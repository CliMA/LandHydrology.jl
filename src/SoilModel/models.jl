export SoilHydrologyModel,
    SoilModel,
    PrescribedTemperatureModel,
    SoilEnergyModel

abstract type AbstractSoilModel end

struct SoilEnergyModel <: AbstractSoilModel
    initial_conditions
    boundary_conditions
end
struct SoilHydrologyModel <: AbstractSoilModel
    initial_conditions
    boundary_conditions
end

struct PrescribedTemperatureModel <: AbstractSoilModel
    profile
end


struct SoilModel{A,B} <: AbstractSoilModel
    energy_model::AbstractSoilModel
    hydrology_model::AbstractSoilModel
    soil_params::A
    earth_params::B
end





