module SoilWaterParameterizations
using Printf
using DocStringExtensions
using UnPack
export vanGenuchten,
    effective_saturation,
    pressure_head,
    hydraulic_conductivity,
    hydrostatic_profile,
    matric_potential,
    viscosity_factor,
    impedance_factor,
    IceImpedance,
    TemperatureDependentViscosity,
    NoEffect,
    AbstractConductivityFactor,
    volumetric_liquid_fraction,
    inverse_matric_potential


"""
    AbstractConductivityFactor{FT <: AbstractFloat}

An abstract type for multiplicative factors that affect the hydraulic conductivity,
increasing or decreasing it relative to the formula that depends only on liquid water content.

Examples include a multiplicative factor due to ice impedance or one due to the
 temperature dependence of the viscosity of water.
"""
abstract type AbstractConductivityFactor{FT <: AbstractFloat} end

"""
    NoEffect{FT} <: AbstractConductivityFactor{FT}

The default type of conductivity factor, used to set a particular multiplicative
factor to unity.
"""
struct NoEffect{FT} <: AbstractConductivityFactor{FT} end

"""
    TemperatureDependentViscosity{FT} <: AbstractConductivityFactor{FT}

The factor type which alters the conductivity to account for the effects of temperature
on the viscosity of water.
"""
Base.@kwdef struct TemperatureDependentViscosity{FT} <:
                   AbstractConductivityFactor{FT}
    "Empirical coefficient"
    γ::FT = FT(2.64e-2)
    "Reference temperature"
    T_ref::FT = FT(288.0)
end

"""
    IceImpedance{FT} <: AbstractConductivityFactor{FT}

The factor type which alters the conductivity to account for the effects of ice
in the pore space.

The empirical factor comes from Lundin (1990).
"""
Base.@kwdef struct IceImpedance{FT} <: AbstractConductivityFactor{FT}
    "Empirical coefficient"
    Ω::FT = FT(7)
end

"""
    impedance_factor(
        imp::NoEffect{FT},
        _...,
    ) where {FT}

Returns a multiplicative factor of 1.0, which does not account for
the ice on conductivity.
"""
function impedance_factor(imp::NoEffect{FT}, _...) where {FT}
    return FT(1.0)
end

"""
    impedance_factor(
        imp::IceImpedance{FT},
        f_i::FT,
    ) where {FT}

Returns the multiplicative factor when an effect due to the fraction of 
ice `f_i` is desired. 
"""
function impedance_factor(imp::IceImpedance{FT}, f_i::FT) where {FT}
    Ω = imp.Ω
    gamma = FT(10.0^(-Ω * f_i))
    return gamma
end

"""
    viscosity_factor(
        vm::NoEffect{FT},
        _...,
    ) where {FT}

Returns a multiplicative factor of 1.0, which does not account for
the temperature dependence of the conductivity.
"""
function viscosity_factor(vm::NoEffect{FT}, _...) where {FT}
    return FT(1.0)
end

"""
    viscosity_factor(
        vm::TemperatureDependentViscosity{FT},
        T::FT,
    ) where {FT}

Returns the multiplicative factor which accounts for
the temperature dependence of the conductivity.
"""
function viscosity_factor(
    vm::TemperatureDependentViscosity{FT},
    T::FT,
) where {FT}
    γ = vm.γ
    T_ref = vm.T_ref
    factor = FT(γ * (T - T_ref))
    Theta = FT(exp(factor))
    return Theta
end

"""
    AbstractsHydraulicsModel{FT <: AbstractFloat}

An abstract type for hydraulics models (soil water retention curve and 
corresponding hydraulic conductivity), such as the van Genuchten model
or the Brooks and Corey model.

The hydraulics model is used in computing the moisture dependence of the hydraulic 
conductivity and in the matric potential. The single hydraulics model 
choice sets both of these.
"""
abstract type AbstractHydraulicsModel{FT <: AbstractFloat} end


"""
    vanGenuchten{FT} <: AbstractHydraulicsModel{FT}


The necessary parameters for the van Genuchten hydraulic model; 
defaults are for loam, setting the residual water fraction to zero.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct vanGenuchten{FT} <: AbstractHydraulicsModel{FT}
    "Pore size distribution index, unitless"
    n::FT
    "Inverse of the air entry potential (1/m)"
    α::FT
    "Exponent parameter - determined as 1-1/n"
    m::FT
    "Residual water content"
    θr::FT
    "Saturated hydraulic conductivity, m/s"
    Ksat::FT
    function vanGenuchten{FT}(;
        n::FT = FT(1.56),
        α::FT = FT(3.6),
        Ksat::FT = FT(2.9e-7),
        θr::FT = FT(0.0),
    ) where {FT}
        new(n, α, FT(1) - FT(1) / FT(n), θr, Ksat)
    end
end

"""
    volumetric_liquid_fraction(
        ϑ_l::FT,
        ν_eff::FT,
    ) where {FT}

Computes the volumetric liquid fraction from the effective porosity and the augmented liquid
fraction.
"""
function volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT) where {FT}
    if ϑ_l < ν_eff
        θ_l = ϑ_l
    else
        θ_l = ν_eff
    end
    return θ_l
end

"""
    matric_potential(hm::vanGenuchten{FT}, S)

Computes the matric potential according the the hydraulics model `hm`,
given an effective saturation `S`.
"""
function matric_potential(hm::vanGenuchten{FT}, S::FT) where {FT}
    @unpack n, α, m = hm
    ψ_m = -((S^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n)
    return ψ_m
end

"""
    effective_saturation(porosity::FT, ϑ_l::FT, θr::FT) 

Computes the effective saturation of the soil, given the augmented liquid fraction
`ϑ_l`, the residual water fraction `θr`, and a porosity.

This porosity can be effective (`ν-θ_i`) or true (ν) depending on the application.
Note that ϑ_l is defined to be larger than `θr`. If `ϑ_l-θr` is negative, 
hydraulic functions that take it as an argument will return 
imaginary numbers, resulting in domain errors.
"""
function effective_saturation(porosity::FT, ϑ_l::FT, θr::FT) where {FT}
    ϑ_l_safe = max(ϑ_l, θr + eps(FT))
    S_l = (ϑ_l_safe - θr) / (porosity - θr)
    return S_l
end

"""
    pressure_head(hm::vanGenuchten{FT}, ϑ_l::FT, ν_eff::FT, S_s::FT)

Computes the pressure head according to the hydraulics model `hm`.

For unsaturated soil (`ϑ_l < ν_eff`, or `S_l_eff < 1`), the matric potential is used. It is
evaluated at an effective saturation computed using the effective porosity. This
ensures that the pressure head is continuous across the saturated and unsaturated 
border. For saturated soil, the pressure head is positive.
"""
function pressure_head(
    hm::vanGenuchten{FT},
    ϑ_l::FT,
    ν_eff::FT,
    S_s::FT,
) where {FT}
    S_l_eff = effective_saturation(ν_eff, ϑ_l, hm.θr)
    if S_l_eff < FT(1.0)
        ψ = matric_potential(hm, S_l_eff)
    else
        ψ = (ϑ_l - ν_eff) / S_s
    end
    return ψ
end

"""
    inverse_matric_potential(
        hm::vanGenuchten{FT},
        ψ::FT
    ) where {FT}

Compute the effective saturation given the matric potential `ψ`, using
the van Genuchten formulation `hm`.
"""
function inverse_matric_potential(hm::vanGenuchten{FT}, ψ::FT) where {FT}
    ψ > 0 && error("Matric potential is positive")
    @unpack n, m, α = hm
    S = (FT(1) + (α * abs(ψ))^n)^(-m)
    return S
end


"""
    hydraulic_conductivity(hm::vanGenuchten{FT}, S::FT, viscosity_f::FT, impedance_f::FT) where {FT}

Computes the hydraulic conductivity given a hydraulics model `hm`, 
the effective saturation, and the additional 
viscosity and ice impedance factors, which must be computed prior 
to evaluating this function.
"""
function hydraulic_conductivity(
    hm::vanGenuchten{FT},
    S::FT,
    viscosity_f::FT,
    impedance_f::FT,
) where {FT}
    @unpack Ksat, m = hm
    if S < FT(1)
        K = sqrt(S) * (FT(1) - (FT(1) - S^(FT(1) / m))^m)^FT(2)
    else
        K = FT(1)
    end
    return K * Ksat * viscosity_f * impedance_f
end

"""
    hydrostatic_profile(hm, z, zmin, ν)

Returns the augmented liquid fraction corresponding to a hydrostatic profile in variably saturated soil,
with the boundary between the saturated and unsaturated zones at `z_∇`.
"""
function hydrostatic_profile(
    hm::vanGenuchten{FT},
    z::FT,
    z_∇::FT,
    ν::FT,
    S_s::FT,
) where {FT}
    @unpack α, m, n, θr = hm
    if z > z_∇
        S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
        ϑ_l = S * (ν - θr) + θr
    else
        ϑ_l = -S_s * (z - z_∇) + ν
    end

    return FT(ϑ_l)
end


end
