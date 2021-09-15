
@testset "soil water params" begin
    FT = Float32
    θr = FT(0.2)
    hm = vanGenuchten{FT}(θr = θr)
    @unpack α, m, n, Ksat = hm
    ν = FT(0.4)
    S_s = FT(1e-2)

    θ = FT.([0.3, 0.4, 0.5])
    S = effective_saturation.(θ; ν = ν, θr = θr)
    @test S ≈ [0.5, 1.0, 1.5]
    @test eltype(S) == FT
    va = -((S[1]^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n)
    ψ = matric_potential.(Ref(hm), S[1:2])
    @test ψ ≈ [va, 0.0]
    @test eltype(ψ) == FT
    p = pressure_head.(Ref(hm), S; ν = ν, S_s = S_s)
    @test p ≈ push!(ψ, FT(10.0))
    @test eltype(p) == FT
    k = hydraulic_conductivity.(Ref(hm), S)
    va = (sqrt(S[1]) * (FT(1) - (FT(1) - S[1]^(FT(1) / m))^m)^FT(2)) * Ksat
    @test k ≈ [va, Ksat, Ksat]
    @test eltype(k) == FT
    z = FT.(Array(-1:0.1:0))
    θ = hydrostatic_profile.(Ref(hm), z, FT(-1.0), ν)
    ψ = matric_potential.(Ref(hm), effective_saturation.(θ; ν = ν, θr = θr))
    h = ψ .+ z
    @test eltype(h) == FT
    @test std(h) < 1e-6
end
