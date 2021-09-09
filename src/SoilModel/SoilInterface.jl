module SoilInterface
import ClimaCore: Fields, Operators, Geometry
using DocStringExtensions
using UnPack

using LandHydrology.Domains: AbstractVerticalDomain, make_function_space
using LandHydrology.Models: AbstractModel
using LandHydrology.BoundaryConditions: AbstractBC, compute_vertical_flux, SoilDomainBC
include("SoilWaterParameterizations.jl")
using .SoilWaterParameterizations
include("SoilHeatParameterizations.jl")
using .SoilHeatParameterizations

import LandHydrology: BoundaryConditions
include("parameters.jl")
include("models.jl")
include("initial_conditions.jl")

end
