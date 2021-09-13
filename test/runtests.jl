using ClimaCore: Fields
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using OrdinaryDiffEq:
    ODEProblem,
    solve,
    CarpenterKennedy2N54,# does not work
    SSPRK33,
    SSPRK73
using LandHydrology
using LandHydrology.Models: make_rhs
using LandHydrology.Domains: Column, make_function_space
using LandHydrology.SoilInterface
using LandHydrology.SoilInterface.SoilWaterParameterizations
using LandHydrology.SoilInterface.SoilHeatParameterizations
using Statistics
using Test
using ArtifactWrappers
using DelimitedFiles

include("test_domains.jl")
include("joint_test.jl")
include("SoilModel/richards_equation.jl")
include("SoilModel/heat_test_interface.jl")
