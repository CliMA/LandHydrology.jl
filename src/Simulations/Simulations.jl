module Simulations

using DiffEqBase

# using ClimaAtmos.Callbacks
using LandHydrology.Models: AbstractModel
using LandHydrology.SoilInterface: make_rhs
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
