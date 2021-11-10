using SurfaceFluxes: surface_conditions, DGScheme
using Thermodynamics: Liquid, q_vap_saturation_generic, PhasePartition, cp_m
using CLIMAParameters.Atmos.Microphysics: D_vapor
using CLIMAParameters.Planet:
    R_v, grav, ρ_cloud_liq, cp_d, R_d, T_0, LH_v0, cp_v
using Printf
export NoBC,
    VerticalFlux,
    Dirichlet,
    FreeDrainage,
    SoilColumnBC,
    SoilComponentBC,
    PrescribedAtmosForcing,
    boundary_fluxes,
    compute_turbulent_surface_fluxes

#### Specific BC types

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
# Fields
$(DocStringExtensions.FIELDS)
"""
struct VerticalFlux{f <: AbstractFloat} <: AbstractBC
    "Scalar flux; positive = aligned with ẑ"
    flux::f
end


"""
    Dirichlet{f} <: AbstractBC

A BC type setting a boundary value of the state `ϑ_l` or
temperature `T`. This can be a function of time.

Note that we apply boundary conditions to temperature even though
volumetric internal energy is the prognostic variable.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Dirichlet{f} <: AbstractBC
    "State Value; (t) -> f(t)"
    state_value::f
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


#### Types for attaching boundary conditions to the soil model

abstract type AbstractFaceBC end

""" 
    SoilComponentBC{ebc <: AbstractBC, hbc <: AbstractBC} <: AbstractFaceBC

A container for holding the boundary conditions for the components of the soil model.

The values must be of type AbstractBC; the two components are energy and hydrology.
Each boundary will have a SoilComponentBC object associated with it.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SoilComponentBC{ebc <: AbstractBC, hbc <: AbstractBC} <:
                   AbstractFaceBC
    "BC for the heat equation"
    energy::ebc = NoBC()
    "BC for ϑ_l"
    hydrology::hbc = NoBC()
end

"""
    PrescribedAtmosForcing{FT <: AbstractFloat} <: AbstractFaceBC

A container for holding the atmosphere state variables that are needed to drive the soil model
via a Monin-Obukhov theory for the surface fluxes.

At present, this boundary condition can only be used with a SoilModel
that includes both volumetric internal energy and volumetric water fractions
as prognostic variables. Sublimation of ice is not currently supported,
nor is precipitation, nor any additional resistance due to a dry layer in
the soil.  Eventually this will support functions of space and time, 
so that reanalysis data can be used to drive the land model. 

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PrescribedAtmosForcing{FT <: AbstractFloat} <: AbstractFaceBC
    "The wind speed parallel to the surface at height z_atm"
    u_atm::FT
    "The potential temperature of the atmosphere at the height z_atm"
    θ_atm::FT
    "The height at which the atmospheric state is measured"
    z_atm::FT
    "The porential temperature scale"
    θ_scale::FT
    "The density of the moist air at the surface"
    ρ_a_sfc::FT
    "The specific humidity of the atmosphere at the height z_atm"
    q_atm::FT
end

"""
    SoilColumnBC{TBC, BBC}

A container holding the Soil BC for each boundary face.

Each field value should be of supertype AbstractFaceBC.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilColumnBC{
    TBC <: Union{SoilComponentBC, PrescribedAtmosForcing},
    BBC <: SoilComponentBC,
}
    "SoilComponentBC for the top of the domain"
    top::TBC
    "Soil BC for the bottom of the domain"
    bottom::BBC
end


function SoilColumnBC(;
    top::Union{SoilComponentBC, PrescribedAtmosForcing} = SoilComponentBC(),
    bottom::SoilComponentBC = SoilComponentBC(),
)
    args = (top, bottom)
    return SoilColumnBC{typeof.(args)...}(args...)
end



###### Methods for expressing any BC in terms of a flux

"""

    interior_values(X::Fields.FieldVector, face::Symbol, cs::Spaces.CenterFiniteDifferenceSpace)::Vector{eltype(X)}

Returns the values of the (center) state variables Y nearest to the boundary
face `face`.
"""
function interior_values(
    X::Fields.FieldVector,
    face::Symbol,
    cs::Spaces.CenterFiniteDifferenceSpace,
)::Vector{eltype(X)}
    @unpack ϑ_l, θ_i, T = X

    N = length(cs.center_local_geometry)
    if face == :top
        return Operators.getidx.([ϑ_l, θ_i, T], Ref(Operators.Interior()), N)
    elseif face == :bottom
        return Operators.getidx.([ϑ_l, θ_i, T], Ref(Operators.Interior()), 1)

    else
        throw(ArgumentError("Expected :top or :bottom"))
    end
end

"""
    boundary_cf_distance(face::Symbol, cs::Spaces.CenterFiniteDifferenceSpace)::Real
Returns the distance from the last center space `cs` node to the nearest face `face`.
"""
function boundary_cf_distance(
    face::Symbol,
    cs::Spaces.CenterFiniteDifferenceSpace,
)::Real
    if face == :top
        return last(cs.face_local_geometry.WJ)
    elseif face == :bottom
        return cs.face_local_geometry.WJ[1]
    else
        throw(ArgumentError("Expected :top or :bottom"))
    end

end

"""
    initialize_boundary_values(X, face::Symbol, model::SoilModel, cs::Spaces.CenterFiniteDifferenceSpace)

Initializes the 2-element arrays (pairs) of boundary values (center, face), `X_cf`,
for each of ϑ_l, θ_i, T.

The initialization sets the boundary face value equal to the center value.
"""
function initialize_boundary_values(
    X,
    face::Symbol,
    model::SoilModel,
    cs::Spaces.CenterFiniteDifferenceSpace,
)
    ϑ_l, θ_i, T = interior_values(X, face, cs)
    θ_l = ϑ_l
    # center, face
    return (ϑ_l = [ϑ_l, ϑ_l], T = [T, T], θ_i = [θ_i, θ_i])
end

"""
    set_boundary_values!(
        X_cf::NamedTuple,
        bc::Dirichlet,
        component::SoilEnergyModel,
        t,
    )

Updates the face element of the boundary value pair for `T` if a Dirichlet
boundary condition on `T` is used.
"""
function set_boundary_values!(
    X_cf::NamedTuple,
    bc::Dirichlet,
    component::SoilEnergyModel,
    t,
)
    X_cf.T[2] = bc.state_value(t)
end

"""
    set_boundary_values!(
        X_cf::NamedTuple,
        bc::Dirichlet,
        component::SoilHydrologyModel,
        t,
    )

Updates the face element of the boundary value pair for `ϑ_l` if a Dirichlet
boundary condition on `ϑ_l` is used.
"""
function set_boundary_values!(
    X_cf::NamedTuple,
    bc::Dirichlet,
    component::SoilHydrologyModel,
    t,
)
    X_cf.ϑ_l[2] = bc.state_value(t)
end


"""
    set_boundary_values!(
        X_cf::NamedTuple,
        bc::AbstractBC,
        component::AbstractSoilComponentModel,
        t,
    )

Does nothing in the case of a non-Dirichlet BC.
"""
function set_boundary_values!(
    X_cf::NamedTuple,
    bc::AbstractBC,
    component::AbstractSoilComponentModel,
    t,
)
    nothing
end

"""
     vertical_flux(bc::VerticalFlux, component::AbstractSoilComponentModel, _...)

Return the vertical flux value at the boundary (flux boundary condition).
"""
function vertical_flux(
    bc::VerticalFlux,
    component::AbstractSoilComponentModel,
    _...,
)
    return bc.flux
end

"""
     vertical_flux(bc::NoBC, _...)

Returns nothing.

Prescribed models, or models that do not require boundary conditions, have `NoBC` as
the boundary type (as opposed to `Dirichlet`, or `VerticalFlux`, e.g.).
 In this case, do not return any flux.
"""
function vertical_flux(bc::NoBC, _...)
    nothing
end

"""
    vertical_flux(
        bc::FreeDrainage,
        component::SoilHydrologyModel,
        X_cf::NamedTuple,
        soil::SoilModel,
        _...,
    )

Returns the vertical flux of water volume at the bottom of the domain
in the case of free drainage.
"""
function vertical_flux(
    bc::FreeDrainage,
    component::SoilHydrologyModel,
    X_cf::NamedTuple,
    soil::SoilModel,
    _...,
)
    @unpack ϑ_l, θ_i, T = X_cf # [center, face]
    ϑ_l = ϑ_l[1]
    θ_i = θ_i[1]
    T = T[1]
    @unpack ν = soil.soil_param_set
    hm = component.hydraulic_model
    @unpack θr = hm

    ν_eff = ν .- θ_i
    θ_l = volumetric_liquid_fraction(ϑ_l, ν_eff)
    T = T[1]
    f_i = θ_i ./ (θ_l .+ θ_i)

    impedance_f = impedance_factor(component.impedance_factor, f_i)
    viscosity_f = viscosity_factor(component.viscosity_factor, T)
    S = effective_saturation(ν, ϑ_l, θr)
    K = hydraulic_conductivity(hm, S, viscosity_f, impedance_f)

    flux = -K # = -K∇h when ∇h = 1. ∇h = 1 -> θ_c = θ_f at the bottom, so use K(θ_c).
    return flux

end

"""
    vertical_flux(
        bc::Dirichlet,
        component::SoilHydrologyModel,
        X_cf::NamedTuple,
        soil::SoilModel,
        dz::AbstractFloat,
        face::Symbol,
    )

Computes the volumetric water flux assuming a Dirichlet condtion
on `ϑ_l` at the boundary.
"""
function vertical_flux(
    bc::Dirichlet,
    component::SoilHydrologyModel,
    X_cf::NamedTuple,
    soil::SoilModel,
    dz::AbstractFloat,
    face::Symbol,
)
    @unpack ϑ_l, θ_i, T = X_cf # [center, face]
    @unpack ν, S_s = soil.soil_param_set
    hm = component.hydraulic_model
    @unpack θr = hm

    ν_eff = ν .- θ_i
    θ_l = volumetric_liquid_fraction.(ϑ_l, ν_eff)

    f_i = θ_i ./ (θ_l .+ θ_i)
    impedance_f = impedance_factor.(Ref(component.impedance_factor), f_i)
    viscosity_f = viscosity_factor.(Ref(component.viscosity_factor), T)
    S = effective_saturation.(ν, ϑ_l, θr)
    K = hydraulic_conductivity.(Ref(hm), S, viscosity_f, impedance_f)

    ψ = pressure_head.(Ref(hm), ϑ_l, ν_eff, S_s)

    flux = -K[2] * (ψ[2] - ψ[1] + dz) / dz
    if face == :bottom # at the bottom
        flux *= -1
    end
    return flux

end

"""
    vertical_flux(
        bc::Dirichlet,
        component::SoilEnergyModel,
        X_cf::NamedTuple,
        soil::SoilModel,
        dz::AbstractFloat,
        face::Symbol,
    )

Computes the energy flux assuming a Dirichlet condtion
on `T` at the boundary.
"""
function vertical_flux(
    bc::Dirichlet,
    component::SoilEnergyModel,
    X_cf::NamedTuple,
    soil::SoilModel,
    dz::AbstractFloat,
    face::Symbol,
)
    @unpack ϑ_l, θ_i, T = X_cf # [center, face]
    @unpack ν, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = soil.soil_param_set
    param_set = soil.earth_param_set
    κ_dry = k_dry(param_set, soil.soil_param_set)

    ν_eff = ν .- θ_i
    θ_l = volumetric_liquid_fraction.(ϑ_l, ν_eff)

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

"""
    boundary_fluxes(
        X,
        bc::SoilComponentBC,
        face::Symbol,
        model::SoilModel,
        cs,
        t,
    )

Returns the boundary flux at the boundary `face` for all
components of the soil model. 

`X` is a FieldVector with components `ϑ_l`, `θ_i`, and `T`. 
The `bc` field contains the boundary condition
types (`Dirichlet`, `VerticalFlux`, etc) for both the energy and hydrology equations.
In the case of explicitly time-varying boundary conditions, we also must supply the current time
`t`.

The remaining fields are required in the case of certain `bc` types, such as Dirichlet.
In this case, we are turning a state boundary condition into a flux boundary condition, and this
requires knowing the side of the domain we are on (encoded in `face`), the values of `X` just interior
to this face, and the center space `cs` itself.
"""
function boundary_fluxes(
    X,
    bc::SoilComponentBC,
    face::Symbol,
    model::SoilModel,
    cs,
    t,
)
    energy = model.energy_model
    hydrology = model.hydrology_model

    X_cf = initialize_boundary_values(X, face, model, cs)
    set_boundary_values!(X_cf, bc.energy, energy, t)
    set_boundary_values!(X_cf, bc.hydrology, hydrology, t)

    dz = boundary_cf_distance(face, cs)
    fρe_int = vertical_flux(bc.energy, energy, X_cf, model, dz, face)
    fϑ_l = vertical_flux(bc.hydrology, hydrology, X_cf, model, dz, face)
    return (fρe_int = fρe_int, fϑ_l = fϑ_l)
end


"""
    function boundary_fluxes(
        X,
        bc::PrescribedAtmosForcing,
        face::Symbol,
        model::SoilModel,
        cs,
        t,
    )

Returns the boundary flux at the boundary `face` for all
components of the soil model in the case of PrescribedAtmosForcing.

`X` is a FieldVector with components `ϑ_l`, `θ_i`, and `T`. 
The `bc` field contains the boundary condition
types (`Dirichlet`, `VerticalFlux`, etc) for both the energy and hydrology equations.
In the case of explicitly time-varying boundary conditions, we also must supply the current time
`t`.

The remaining fields are required in the case of certain `bc` types, such as Dirichlet.
In this case, we are turning a state boundary condition into a flux boundary condition, and this
requires knowing the side of the domain we are on (encoded in `face`), the values of `X` just interior
to this face, and the center space `cs` itself.
"""
function boundary_fluxes(
    X,
    bc::PrescribedAtmosForcing,
    face::Symbol,
    model::SoilModel,
    cs,
    t,
)
    if face != :top
        error(
            "Prescribed atmosphere-driven boundary conditions are only valid at the top of the soil column.",
        )
    end

    energy = model.energy_model
    hydrology = model.hydrology_model
    ϑ_l, θ_i, T = interior_values(X, face, cs)
    fρe_int, fϑ_l =
        compute_turbulent_surface_fluxes(energy, hydrology, model, ϑ_l, θ_i, T)
    return (fρe_int = fρe_int, fϑ_l = fϑ_l)
end

"""
    compute_turbulent_surface_fluxes(
        energy::SoilEnergyModel,
        hydrology::SoilHydrologyModel,
        model::SoilModel,
        ϑ_l::FT,
        θ_i::FT,
        T::FT,
    ) where {FT}

Returns the surface fluxes resulting from the conditions
at the surface (obtained from the soil model) and from the
prescribed atmosphere variables.

This function employs Monin-Obukhov similarity theory. For more
details, see github.com/Clima/SurfaceFluxes.jl.
"""
function compute_turbulent_surface_fluxes(
    energy::SoilEnergyModel,
    hydrology::SoilHydrologyModel,
    model::SoilModel,
    ϑ_l::FT,
    θ_i::FT,
    T::FT,
) where {FT}

    # Atmospheric state at height z_atm, and surface moist air density
    atmos_conditions = model.boundary_conditions.top
    @unpack θ_scale, z_atm, ρ_a_sfc, u_atm, q_atm, θ_atm = atmos_conditions
    x_in = [u_atm, θ_atm, q_atm]

    # Parameters, including surface roughness z_0
    @unpack ν, z_0m, z_0s = model.soil_param_set
    z_0 = [z_0m, z_0s, z_0s]
    param_set = model.earth_param_set

    # Compute specific humidity in the soil pore air near the surface
    q_sat = q_vap_saturation_generic(param_set, T, ρ_a_sfc, Liquid())
    Rv = R_v(param_set)
    g = grav(param_set)
    hm = hydrology.hydraulic_model
    ν_eff = ν - θ_i
    θ_l = volumetric_liquid_fraction(ϑ_l, ν_eff)
    S_l_eff = min(effective_saturation(ν_eff, θ_l, hm.θr), 1.0)
    ψ = matric_potential(hm, S_l_eff)
    correction = exp(g * ψ / Rv / T)
    q_surf = q_sat * correction
    # The surface state
    x_s = [FT(0.0), T, q_surf]

    Lmo_guess = [
        FT(100.0) * z_atm,
        atmos_conditions.u_atm,
        atmos_conditions.θ_atm,
        q_atm,
    ] # save from, previous step?

    conditions = surface_conditions(
        param_set,
        Lmo_guess,
        x_in,
        x_s,
        z_0,
        θ_scale,
        z_atm,
        DGScheme(),
    )

    (ustar, tstar, qstar) = conditions.x_star
    qp = PhasePartition(q_surf) # assumes all moisture is in vapor phase
    cpm = cp_m(param_set, qp)
    _T_ref = FT(T_0(param_set))
    h_d = cp_d(param_set) * (T - _T_ref) + R_d(param_set) * _T_ref

    # Fluxes of energy and water volume
    E = -ρ_a_sfc * ustar * qstar
    dry_static_energy_flux = -cpm * ρ_a_sfc * ustar * tstar - h_d * E # ignoring change in height
    vapor_static_energy_flux =
        FT(cp_v(param_set) * (T - _T_ref) + LH_v0(param_set)) * E# ignoring change in height
    Ẽ = E / ρ_cloud_liq(param_set) # the soil model needs a volume flux
    heat_flux = dry_static_energy_flux + vapor_static_energy_flux
    return heat_flux, Ẽ
end
