module SoilWaterParameterizations
using Printf
using DocStringExtensions
using UnPack
export vanGenuchten,
    effective_saturation,
    pressure_head,
    hydraulic_conductivity,
    hydrostatic_profile,
    matric_potential


"""
    AbstractsHydraulicsModel{FT <: AbstractFloat}

Hydraulics model is used in the moisture factor in hydraulic 
conductivity and in the matric potential. The single hydraulics model 
choice sets both of these.
"""
abstract type AbstractHydraulicsModel{FT <: AbstractFloat} end


"""
    vanGenuchten{FT} <: AbstractHydraulicsModel{FT}


The necessary parameters for the van Genuchten hydraulic model; 
defaults are for Yolo light clay.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct vanGenuchten{FT} <: AbstractHydraulicsModel{FT}
    "Exponent parameter - used in matric potential"
    n::FT
    "used in matric potential. The inverse of this carries units in 
     the expression for matric potential (specify in inverse meters)."
    α::FT
    "Exponent parameter - determined by n, used in hydraulic conductivity"
    m::FT
    "Residual water content"
    θr::FT
    "Saturated hydraulic conductivity, m/s"
    Ksat::FT
    function vanGenuchten{FT}(;
        n::FT = FT(1.43),
        α::FT = FT(2.6),
        Ksat::FT = FT(1e-6),
        θr = FT(0.0),
    ) where {FT}
        new(n, α, FT(1) - FT(1) / FT(n), θr, Ksat)
    end
end



"""
    function matric_potential(hm::vanGenuchten{FT}, S)
"""
function matric_potential(hm::vanGenuchten{FT}, S::FT) where {FT}
    @unpack n, α, m = hm
    ψ_m = -((S^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n)
    return ψ_m
end

"""
    function effective_saturation(θ::FT; ν::FT = ν, θr::FT = θr)
"""
function effective_saturation(θ::FT; ν::FT = ν, θr::FT = θr) where {FT}
    if θ < θr
        println("Effective saturation is negative")
        println(θ)
    end
    return (θ - θr) / (ν - θr)
end

"""
    function pressure_head(hm::vanGenuchten{FT}, S::FT; ν::FT = ν, S_s::FT = S_s)
"""
function pressure_head(
    hm::vanGenuchten{FT},
    S::FT;
    ν::FT = ν,
    S_s::FT = S_s,
) where {FT}
    @unpack θr = hm
    if S < FT(0)
        println("Effective saturation is negative")
        println(S)
    end
    if S < FT(1)
        ψ = matric_potential(hm, S)
    else
        θ = S * (ν - θr) + θr
        ψ = (θ - ν) / S_s
    end

    return ψ
end

"""
    function hydraulic_conductivity(hm::vanGenuchten{FT}, S::FT)
"""
function hydraulic_conductivity(hm::vanGenuchten{FT}, S::FT) where {FT}
    @unpack Ksat, m = hm
    if S < FT(1)
        K = sqrt(S) * (FT(1) - (FT(1) - S^(FT(1) / m))^m)^FT(2)
    else
        K = FT(1)
    end
    return K * Ksat
end

"""
    function hydrostatic_profile(hm, z, zmin, ν)
"""
function hydrostatic_profile(
    hm::vanGenuchten{FT},
    z::FT,
    zmin::FT,
    ν::FT,
) where {FT}
    @unpack α, m, n, θr = hm
    S = FT((FT(1) + (α * (z - zmin))^n)^(-m))
    return FT(S * (ν - θr) + θr)
end


end
