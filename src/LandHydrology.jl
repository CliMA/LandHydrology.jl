module LandHydrology
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end

include("Domains/Domains.jl")
include("Models.jl")
include(joinpath("SoilModel", "SoilInterface.jl"))
include("Simulations/Simulations.jl")

end # module
