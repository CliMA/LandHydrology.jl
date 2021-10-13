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
# Fields
$(DocStringExtensions.FIELDS)
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
    SoilComponentBC{ebc <: AbstractBC, hbc <: AbstractBC}

A container for holding the boundary conditions for the components of the soil model.

The values must be of type AbstractBC; the two components are energy and hydrology.
Each boundary will have a SoilComponentBC object associated with it.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SoilComponentBC{ebc <: AbstractBC, hbc <: AbstractBC}
    "BC for the heat equation"
    energy::ebc = NoBC()
    "BC for ϑ_l"
    hydrology::hbc = NoBC()
end


"""
    SoilDomainBC{D, TBC, BBC}

A container holding the SoilComponentBC for each boundary face.

Each field value should be of type SoilComponentBC. This doesn't do 
what we want. Ideally the fields would change depending on the domain - 
e.g. for a Column, they are top and bottom, but for a 3D domain, they might be
top, bottom, xleft, xright, yleft, yright.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilDomainBC{D, TBC, BBC}
    "SoilComponentBC for the top of the domain"
    top::TBC
    "SoilComponentBC for the bottom of the domain"
    bottom::BBC
end


function SoilDomainBC(
    ::Column{FT};
    top::SoilComponentBC = SoilComponentBC(),
    bottom::SoilComponentBC = SoilComponentBC(),
) where {FT}
    args = (top, bottom)
    return SoilDomainBC{Column{FT}, typeof.(args)...}(args...)
end


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
