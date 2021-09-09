### Soil specific boundary condition methods and types

"""
    FreeDrainage <: AbstractBC

A BC type for use with the SoilHydrologyModel (Richards Equation),
at the bottom of the domain, setting 

``
    ∇h = 1.
``

"""
struct FreeDrainage <: AbstractBC end


function BoundaryConditions.compute_vertical_flux(bc::FreeDrainage, Y)# this only make sense to use at the bottom, but the user should know this. 
    (Y_hydro, Y_energy) = Y.x
    @unpack ϑ_l, θ_i = Y_hydro
    θ_l = ϑ_l
    S = effective_saturation.(θ_l; ν = ν, θr = θr)
    K = hydraulic_conductivity.(S; vgm = vgm, ksat = ksat)
    flux = -parent(K)[1] # = -K∇h when ∇h = 1. ∇h = 1 -> θ_c = θ_f at the bottom, so use K(θ_c).
    return Geometry.Cartesian3Vector(flux)
end
