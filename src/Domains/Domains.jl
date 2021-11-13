module Domains

using ClimaCore
using Printf
using DocStringExtensions

import ClimaCore: Meshes, Spaces, Topologies, Geometry, Fields

"""
    AbstractDomain

An abstract type for domains.
"""
abstract type AbstractDomain end

"""
    AbstractVerticalDomain{FT}

An abstract type for vertical domains, using the floating point
type FT. 
"""
abstract type AbstractVerticalDomain{FT} <: AbstractDomain end

include("domain.jl")

export AbstractDomain
export AbstractVerticalDomain
export Column
export make_function_space

end # module
