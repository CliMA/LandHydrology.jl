module LandHydrology

include(joinpath("Domains", "Domains.jl"))
include("Models.jl")
include(joinpath("SoilModel", "SoilInterface.jl"))
include(joinpath("Simulations", "Simulations.jl"))

end # module
