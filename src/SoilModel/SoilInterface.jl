module SoilInterface
import ClimaCore: Fields, Operators, Geometry
using DocStringExtensions
using UnPack

using LandHydrology.Domains: AbstractVerticalDomain, make_function_space, Column
using LandHydrology.Models: AbstractModel
import LandHydrology: Models, LandHydrology, set_initial_state

include("SoilWaterParameterizations.jl")
using .SoilWaterParameterizations
include("SoilHeatParameterizations.jl")
using .SoilHeatParameterizations

include("parameters.jl")
include("models.jl")
include("initial_conditions.jl")

end
