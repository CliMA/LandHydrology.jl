export VerticalFlux, SoilDomainBC, SoilComponentBC, NoBC, compute_vertical_flux
abstract type AbstractBC end

"""
    NoBC <: AbstractBC

The BC type to be used when equations do not require boundary conditions,
e.g. for Prescribed Models.
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

"""
    FreeDrainage <: AbstractBC

A BC type for use with the SoilHydrologyModel (Richards Equation),
at the bottom of the domain, setting 

``
    ∇h = 1.
``

"""
struct FreeDrainage <: AbstractBC end

Base.@kwdef struct SoilDomainBC{TBC, BBC}
    top::TBC = SoilComponentBC()
    bottom::BBC = SoilComponentBC()
end

Base.@kwdef struct SoilComponentBC{ebc <: AbstractBC, hbc <: AbstractBC}
    energy::ebc = NoBC()
    hydrology::hbc = NoBC()
end


function compute_vertical_flux(bc::VerticalFlux, _...)
    return Geometry.Cartesian3Vector(bc.vertical_flux)
end
#unclear if we need this
function compute_vertical_flux(bc::NoBC, _...)
    return nothing
end

function compute_vertical_flux(bc::FreeDrainage, Y)# this only make sense to use at the bottom, but the user should know this. 
    (Y_hydro, Y_energy) = Y.x
    @unpack ϑ_l, θ_i = Y_hydro
    θ_l = ϑ_l
    S = effective_saturation.(θ_l; ν = ν, θr = θr)
    K = hydraulic_conductivity.(S; vgm = vgm, ksat = ksat)
    flux = -parent(K)[1] # = -K∇h when ∇h = 1. ∇h = 1 -> θ_c = θ_f at the bottom, so use K(θ_c).
    return Geometry.Cartesian3Vector(flux)
end
