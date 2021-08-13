module SoilWaterParameterizations
using Printf

export effective_saturation, pressure_head, hydraulic_conductivity, hydrostatic_profile, matric_potential


function matric_potential(S; vgn = vgn, vgα = vgα, vgm = vgm)
    FT = eltype(S)
    ψ_m = -((S^(-FT(1) / vgm) - FT(1)) * vgα^(-vgn))^(FT(1) / vgn)
    return ψ_m
end


function effective_saturation(θ; ν = ν, θr = θr)
    
    if θ < θr
        println("Effective saturation is negative")
        println(θ)
    end
    return (θ-θr)/(ν-θr)
end

function pressure_head(S; vgn = vgn, vgα = vgα, vgm = vgm, ν = ν, θr = θr, S_s = S_s)
    FT = eltype(S)
    if S < FT(0)
        println("Effective saturation is negative")
        println(S)
    end
    if S < FT(1)
            ψ = matric_potential(S;vgn = vgn, vgα = vgα, vgm = vgm)
    else
        θ = S* (ν-θr) + θr 
        ψ = (θ-ν)/S_s
    end
    
    return ψ
end


function hydraulic_conductivity(S; vgm = vgm, ksat = ksat)
    FT = eltype(S)
    if S < FT(1)
        K = sqrt(S) * (FT(1) - (FT(1) - S^(FT(1) / vgm))^vgm)^FT(2)
    else
        K = FT(1)
    end
        return K*ksat
end
function hydrostatic_profile(z, zmin, porosity, n, α, θr)
    FT = eltype(z)
    m = FT(1 - 1 / n)
    S = FT((FT(1) + (α * (z - zmin))^n)^(-m))
    return FT(S * (porosity-θr)+θr)
end


end
