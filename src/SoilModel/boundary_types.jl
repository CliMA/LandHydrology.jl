export VerticalFlux,
    Dirichlet, FreeDrainage, SoilColumnBC, SoilComponentBC, NoBC
abstract type AbstractBC end
abstract type AbstractSoilDomainBC end
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
    SoilColumnBC{TBC, BBC} <: AbstractSoilDomainBC

A container holding the SoilComponentBC for each boundary face.

Each field value should be of type SoilComponentBC.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilColumnBC{TBC, BBC} <: AbstractSoilDomainBC
    "SoilComponentBC for the top of the domain"
    top::TBC
    "SoilComponentBC for the bottom of the domain"
    bottom::BBC
end


function SoilColumnBC(;
    top::SoilComponentBC = SoilComponentBC(),
    bottom::SoilComponentBC = SoilComponentBC(),
)
    args = (top, bottom)
    return SoilColumnBC{typeof.(args)...}(args...)
end

