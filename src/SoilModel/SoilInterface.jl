module SoilInterface
import ClimaCore: Fields, Operators, Geometry, Spaces
using CLIMAParameters
using CLIMAParameters.Planet: T_0

using DocStringExtensions
using UnPack

using LandHydrology: EarthParameterSet, LandHydrologyModel
using LandHydrology.Domains: AbstractVerticalDomain, make_function_space, Column
using LandHydrology.Models: AbstractModel
import LandHydrology: Models, make_tendency_terms, make_update_aux, initialize_states

include("SoilWaterParameterizations.jl")
using .SoilWaterParameterizations
include("SoilHeatParameterizations.jl")
using .SoilHeatParameterizations

include("parameters.jl")
include("models.jl")
include("initial_conditions.jl")

end
