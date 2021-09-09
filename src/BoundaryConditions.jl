module BoundaryConditions
export NoBC, VerticalFlux, compute_vertical_flux,SoilDomainBC, SoilComponentBC, AbstractBC
abstract type AbstractBC end

"""
    NoBC <: AbstractBC

The BC type to be used when equations do not require boundary conditions,
e.g. for Prescribed Models, or ODEs.
"""
struct NoBC <: AbstractBC end

"""
    VerticalFlux{f <: AbstractFloat} <: AbstractBC

The BC type to be used for prescribed vertical boundary fluxes. The flux is
assumed to be of the form

``
F = f ẑ
``

where f is the value supplied by the user (currently a constant).
"""
struct VerticalFlux{f <: AbstractFloat} <: AbstractBC
    "Scalar flux; positive = aligned with ẑ"
    vertical_flux::f
end

Base.@kwdef struct SoilDomainBC{TBC, BBC}
    top::TBC = SoilComponentBC()
    bottom::BBC = SoilComponentBC()
end

Base.@kwdef struct SoilComponentBC{ebc <: AbstractBC, hbc <: AbstractBC}
    energy::ebc = NoBC()
    hydrology::hbc = NoBC()
end

function compute_vertical_flux(bc::AbstractBC, _...)
end
function compute_vertical_flux(bc::NoBC, _...)
    return nothing
end
function compute_vertical_flux(bc::VerticalFlux, _...)
    return Geometry.Cartesian3Vector(bc.vertical_flux)
end
end
