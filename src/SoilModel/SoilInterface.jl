module SoilInterface
import ClimaCore: Fields, Operators, Geometry
using DocStringExtensions
using UnPack

using LandHydrology.Domains: AbstractVerticalDomain, make_function_space
using LandHydrology.Models: AbstractModel
include("SoilWaterParameterizations.jl")
using .SoilWaterParameterizations
include("SoilHeatParameterizations.jl")
using .SoilHeatParameterizations

include("parameters.jl")
include("models.jl")
include("initial_conditions.jl")

end
