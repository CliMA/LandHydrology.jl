module Simulations

using DiffEqBase

# using ClimaAtmos.Callbacks
using LandHydrology: AbstractLandModel,make_rhs, default_land_initial_conditions
using ClimaCore: Fields
using DocStringExtensions

import DiffEqBase: step!

"""
    AbstractSimulation
"""
abstract type AbstractSimulation end

include("simulation.jl")

export Simulation
export step!
export run!

end # module
