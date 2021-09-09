module Domains

using ClimaCore
using Printf

import ClimaCore: Meshes, Spaces, Topologies

"""
    AbstractDomain
"""
abstract type AbstractDomain end

"""
    AbstractVerticalDomain
"""
abstract type AbstractVerticalDomain{FT} <: AbstractDomain end

include("domain.jl")

export AbstractDomain
export AbstractVerticalDomain
export Column
export make_function_space

end # module
