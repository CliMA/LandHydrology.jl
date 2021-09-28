
@testset "soil water params" begin
    FT = Float32
    θr = FT(0.2)
    hm = vanGenuchten{FT}(θr = θr)
    @unpack α, m, n, Ksat = hm
    ν = FT(0.4)
    S_s = FT(1e-2)

    # Effective saturation
    θ = FT.([0.3, 0.4, 0.5])
    S = effective_saturation.(ν, θ, θr)
    @test S ≈ [0.5, 1.0, 1.5]
    @test eltype(S) == FT
    @test_throws Exception effective_saturation(0.5, -1.0, 0.0)

    # Matric Potential and inverse
    va = -((S[1]^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n)
    ψ = matric_potential.(Ref(hm), S[1:2])
    @test inverse_matric_potential.(Ref(hm), ψ) ≈ S[1:2]
    @test_throws Exception inverse_matric_potential(hm, 1)
    @test ψ ≈ [va, 0.0]
    @test eltype(ψ) == FT

    #Pressure head
    p = pressure_head.(Ref(hm), θ, ν, S_s)
    @test p ≈ push!(ψ, FT(10.0))
    @test eltype(p) == FT

    # Hydraulic K

    CF = NoEffect{FT}()
    vf = viscosity_factor(CF)
    impf = impedance_factor(CF)
    k = hydraulic_conductivity.(Ref(hm), S, vf, impf)
    va = (sqrt(S[1]) * (FT(1) - (FT(1) - S[1]^(FT(1) / m))^m)^FT(2)) * Ksat
    @test k ≈ [va, Ksat, Ksat]
    @test eltype(k) == FT


    # Impedance factor
    ImpF = IceImpedance{FT}()
    @test impedance_factor(ImpF, FT(1.0)) ≈ 1e-7

    # Viscosity Factor
    VisF = TemperatureDependentViscosity{FT}()
    T = FT.([278.0, 288.0, 298.0])
    @test viscosity_factor.(Ref(VisF), T) ≈ exp.(VisF.γ .* (T .- VisF.T_ref))

    # Hydrostatic Profile
    z = FT.(Array(-1:0.1:0))
    θ = hydrostatic_profile.(Ref(hm), z, FT(-0.5), ν, S_s)
    ψ = pressure_head.(Ref(hm), θ, ν, S_s)
    h = ψ .+ z
    @test eltype(h) == FT
    @test std(h) < 1e-6

    # Volumetric Liquid Fraction
    vlf = volumetric_liquid_fraction.(FT.([0.25, 0.5, 0.75]), FT(0.5))
    @test vlf ≈ FT.([0.25, 0.5, 0.5])




end
