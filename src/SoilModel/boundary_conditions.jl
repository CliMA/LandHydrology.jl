export VerticalFlux,
    Dirichlet, FreeDrainage, SoilDomainBC, SoilComponentBC, NoBC
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
    flux::f
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

struct Dirichlet{f} <: AbstractBC
    "State Value; (t) -> f(t)"
    state_value::f
end


Base.@kwdef struct SoilComponentBC{ebc <: AbstractBC, hbc <: AbstractBC}
    energy::ebc = NoBC()
    hydrology::hbc = NoBC()
end

struct SoilDomainBC{D, TBC, BBC}
    top::TBC
    bottom::BBC
end


## It is necessary that the fields of DomainBC match the x3 boundary symbols in Column object
function SoilDomainBC(
    ::Column{FT};
    top::SoilComponentBC = SoilComponentBC(),
    bottom::SoilComponentBC = SoilComponentBC(),
) where {FT}
    args = (top, bottom)
    return SoilDomainBC{Column{FT}, typeof.(args)...}(args...)
end



function get_values(Y, face::Symbol)
    @unpack ϑ_l, θ_i, ρe_int = Y
    if face == :top
        return last(parent(ϑ_l)), last(parent(θ_i)), last(parent(ρe_int))
    elseif face == :bottom
        return first(parent(ϑ_l)), first(parent(θ_i)), first(parent(ρe_int))
    else
        throw(ArgumentError("Expected :top or :bottom"))
    end
end

function get_distance(face::Symbol, cs)
    FT = eltype(cs.mesh)
    if face == :top
        return last(cs.face_local_geometry.WJ)
    elseif face == :bottom
        return cs.face_local_geometry.WJ[1]
    else
        throw(ArgumentError("Expected :top or :bottom"))
    end

end


function initialize_boundary_values(Y, face::Symbol, model::SoilModel)
    ϑ_l, θ_i, ρe_int = get_values(Y, face)
    θ_l = ϑ_l
    @unpack ρc_ds = model.soil_param_set
    param_set = model.earth_param_set
    ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρc_ds, param_set)
    T = temperature_from_ρe_int(ρe_int, θ_i, ρc_s, param_set)
    # center, face
    return (ϑ_l = [ϑ_l, ϑ_l], T = [T, T], θ_i = [θ_i, θ_i])
end


function set_boundary_values!(
    Y_cf::NamedTuple,
    bc::Dirichlet,
    component::SoilEnergyModel,
    t,
)
    Y_cf.T[2] = bc.state_value(t)
end
function set_boundary_values!(
    Y_cf::NamedTuple,
    bc::Dirichlet,
    component::SoilHydrologyModel,
    t,
)
    Y_cf.ϑ_l[2] = bc.state_value(t)
end
function set_boundary_values!(
    Y_cf::NamedTuple,
    bc::AbstractBC,
    component::AbstractSoilComponentModel,
    t,
)
    nothing
end

function get_vertical_flux(bc::VerticalFlux, component::SoilEnergyModel, _...)
    return bc.flux
end
function get_vertical_flux(
    bc::VerticalFlux,
    component::SoilHydrologyModel,
    _...,
)
    return bc.flux
end
function get_vertical_flux(bc::NoBC, _...) # PrescribedModels have NoBC() type
    nothing
end


function get_vertical_flux(
    bc::FreeDrainage,
    component::SoilHydrologyModel,
    Y_cf::NamedTuple,
    soil::SoilModel,
    _...,
)
    @unpack ϑ_l, θ_i, T = Y_cf # [center, face]
    θ_l = ϑ_l[1]
    θ_i = θ_i[1]
    T = T[1]
    @unpack ν, vgm, ksat, θr = soil.soil_param_set
    S = effective_saturation(θ_l; ν = ν, θr = θr)
    K = hydraulic_conductivity(S; vgm = vgm, ksat = ksat)
    flux = -K # = -K∇h when ∇h = 1. ∇h = 1 -> θ_c = θ_f at the bottom, so use K(θ_c).
    return flux

end

function get_vertical_flux(
    bc::Dirichlet,
    component::SoilHydrologyModel,
    Y_cf::NamedTuple,
    soil::SoilModel,
    dz::AbstractFloat,
    face::Symbol,
)
    @unpack ϑ_l, θ_i, T = Y_cf # [center, face]
    θ_l = ϑ_l

    @unpack ν, vgα, vgn, vgm, ksat, θr, S_s = soil.soil_param_set
    S = effective_saturation.(θ_l; ν = ν, θr = θr)
    K = hydraulic_conductivity.(S; vgm = vgm, ksat = ksat)
    ψ =
        pressure_head.(
            S;
            vgn = vgn,
            vgα = vgα,
            vgm = vgm,
            ν = ν,
            θr = θr,
            S_s = S_s,
        )
    flux = -K[2] * (ψ[2] - ψ[1] + dz) / dz
    if face == :bottom # at the bottom
        flux *= -1
    end
    return flux

end

function get_vertical_flux(
    bc::Dirichlet,
    component::SoilEnergyModel,
    Y_cf::NamedTuple,
    soil::SoilModel,
    dz::AbstractFloat,
    face::Symbol,
)
    @unpack ϑ_l, θ_i, T = Y_cf # [center, face]
    @unpack ν, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = soil.soil_param_set
    param_set = soil.earth_param_set
    κ_dry = k_dry(param_set, soil.soil_param_set)

    θ_l = ϑ_l
    S_r = relative_saturation.(θ_l, θ_i, ν)
    kersten = kersten_number.(θ_i, S_r, Ref(soil.soil_param_set))
    κ_sat =
        saturated_thermal_conductivity.(θ_l, θ_i, κ_sat_unfrozen, κ_sat_frozen)
    κ = thermal_conductivity.(κ_dry, kersten, κ_sat) # at face

    flux = -κ[2] * (T[2] - T[1]) / dz
    if face == :bottom # at the bottom
        flux *= -1
    end

    return flux
end

# This last function is the only one needed outside of this. Called inside
# the RHS function to figure out boundary fluxes

function return_fluxes(
    Y,
    bc::SoilComponentBC,
    face::Symbol,
    model::SoilModel,
    cs,
    t,
)
    energy = model.energy_model
    hydrology = model.hydrology_model

    Y_cf = initialize_boundary_values(Y, face, model)
    set_boundary_values!(Y_cf, bc.energy, energy, t)
    set_boundary_values!(Y_cf, bc.hydrology, hydrology, t)

    dz = get_distance(face, cs)
    fρe_int = get_vertical_flux(bc.energy, energy, Y_cf, model, dz, face)
    fϑ_l = get_vertical_flux(bc.hydrology, hydrology, Y_cf, model, dz, face)
    return (fρe_int = fρe_int, fϑ_l = fϑ_l)
end
