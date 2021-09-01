export VerticalFlux,SoilDomainBC, SoilComponentBC, NoBC, compute_vertical_flux
abstract type AbstractBC end

struct NoBC <: AbstractBC end

struct VerticalFlux{f <: AbstractFloat} <: AbstractBC
    "Scalar flux; positive = aligned with ẑ"
    vertical_flux::f
end

struct FreeDrainage <: AbstractBC
end

Base.@kwdef struct SoilDomainBC{TBC, BBC}
    top::TBC = SoilComponentBC()
    bottom::BBC = SoilComponentBC()
end

Base.@kwdef struct SoilComponentBC{ebc<: AbstractBC, hbc<:AbstractBC}
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
