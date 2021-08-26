module SoilInterface
import ClimaCore:    Fields, Operators
using DocStringExtensions
using UnPack

include("SoilWaterParameterizations.jl")
using .SoilWaterParameterizations
include("SoilHeatParameterizations.jl")
using .SoilHeatParameterizations

include("models.jl")
include("initial_conditions.jl")

end
