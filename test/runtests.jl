using ClimaCore: Fields
using UnPack
using OrdinaryDiffEq:
    ODEProblem,
    solve,
    CarpenterKennedy2N54,# does not work
    SSPRK33,
    SSPRK73
using LandHydrology
using LandHydrology: EarthParameterSet
const param_set = EarthParameterSet()
using LandHydrology.Domains: Column, make_function_space
using LandHydrology.SoilInterface
using LandHydrology.SoilInterface.SoilWaterParameterizations
using LandHydrology.SoilInterface.SoilHeatParameterizations
using Statistics
using Test
using ArtifactWrappers
using DelimitedFiles

include("test_domains.jl")
@testset "Soil Model" begin
    @info "Testing LandHydrology Soil model"
    include("SoilModel/coupled.jl")
    include("SoilModel/richards_equation.jl")
    include("SoilModel/heat_test_interface.jl")
    include("SoilModel/test_rhs.jl")
    include("SoilModel/test_water_parameterizations.jl")
    include("SoilModel/test_heat_parameterizations.jl")
end
