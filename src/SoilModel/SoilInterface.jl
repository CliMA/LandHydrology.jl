module SoilInterface
import ClimaCore: Fields, Operators, Geometry, Spaces
using CLIMAParameters
using CLIMAParameters.Planet: ρ_cloud_liq, ρ_cloud_ice, T_freeze, grav, LH_f0
using DocStringExtensions
using UnPack

using LandHydrology.Domains: AbstractVerticalDomain, make_function_space, Column
using LandHydrology.Models: AbstractModel, AbstractLandSource
include("SoilWaterParameterizations.jl")
using .SoilWaterParameterizations
include("SoilHeatParameterizations.jl")
using .SoilHeatParameterizations

include("parameters.jl")
include("models.jl")
include("sources.jl")
include("boundary_conditions.jl")
include("right_hand_side.jl")
include("initial_conditions.jl")

end
