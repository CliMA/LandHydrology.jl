using ClimaCore: Fields
using CLIMAParameters
using CLIMAParameters.Planet: T_0
using UnPack
using OrdinaryDiffEq:
    ODEProblem,
    solve,
    SSPRK33,
    SSPRK73
using LandHydrology
using LandHydrology: make_rhs, LandHydrologyModel
using LandHydrology.Models: default_initial_conditions, make_tendency_terms, make_update_aux, initialize_states, NotIncluded
using LandHydrology: EarthParameterSet
const param_set = EarthParameterSet()
using LandHydrology.Domains: Column, make_function_space
using LandHydrology.SoilInterface
using LandHydrology.SoilInterface.SoilWaterParameterizations
using LandHydrology.SoilInterface.SoilHeatParameterizations
using LandHydrology.Simulations
using LandHydrology.SurfaceWater
using Statistics
using Test
using ArtifactWrappers
using DelimitedFiles

include("test_domains.jl")

@testset "LandHydrology unit tests" begin
end

@testset "Soil+Surface" begin
    @info "Testing the soil and surface model together"
    include("SurfaceWater/joint.jl")
end

@testset "Soil Model" begin
    @info "Testing LandHydrology Soil model"
    include("SoilModel/coupled.jl")
    include("SoilModel/richards_equation.jl")
    include("SoilModel/heat_test_interface.jl")
    include("SoilModel/test_rhs.jl")
    include("SoilModel/test_water_parameterizations.jl")
    include("SoilModel/test_heat_parameterizations.jl")
end
