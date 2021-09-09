import ClimaCore: Fields, Domains, Meshes, Operators, Spaces, Geometry
using OrdinaryDiffEq:
    ODEProblem,
    solve,
    SSPRK33,
    Rosenbrock23,
    Tsit5,
    SSPRK432,
    Feagin14,
    TsitPap8,
    CarpenterKennedy2N54
using DelimitedFiles
using UnPack
using LandHydrology
using LandHydrology.SoilInterface: SoilParams
using LandHydrology.SoilInterface.SoilWaterParameterizations
using ArtifactWrappers
using Test

const FT = Float64


abstract type bc end

struct θDirichlet{FT <: AbstractFloat} <: bc
    θvalue::FT
end
struct FreeDrainage <: bc end



function compute_richards_rhs!(
    dθl,
    θl,
    p,
    top::θDirichlet,
    bottom::FreeDrainage,
)
    # θ_top = something, K∇h_bottom = K_bottom (∇h = 1). If ∇h = 1, θ_b = θ_f, so K = Kbottom at the face.

    sp = p[2]
    zc = p[1]
    @unpack ν, vgn, vgα, vgm, ksat, θr = sp
    S = effective_saturation.(θl; ν = ν, θr = θr)
    K = hydraulic_conductivity.(S; vgm = vgm, ksat = ksat)
    ψ = matric_potential.(S; vgn = vgn, vgα = vgα, vgm = vgm)
    h = ψ .+ zc


    θ_top = top.θvalue
    S_top = effective_saturation.(θ_top; ν = ν, θr = θr)
    h_top = matric_potential(S_top; vgn = vgn, vgα = vgα, vgm = vgm)#ztop = 0
    K_top = hydraulic_conductivity(S_top; vgm = vgm, ksat = ksat)
    bc_ktop = Operators.SetValue(K_top)
    If = Operators.InterpolateC2F(; top = bc_ktop)

    bc_t = Operators.SetValue(h_top)
    gradc2f = Operators.GradientC2F(top = bc_t) # set value on h_top

    bc_b = Operators.SetValue(Geometry.Cartesian3Vector(FT(1) * parent(K)[1]))
    divf2c = Operators.DivergenceF2C(bottom = bc_b) # set value on K∇h at bottom



    return @. dθl = divf2c(If(K) * gradc2f(h))
end

@testset "Richards sand 1" begin
    n = 150
    z₀ = FT(-1.5)
    z₁ = FT(0)
    ksat = FT(34 / (3600 * 100))
    vgn = FT(3.96)
    vgα = FT(2.7)
    vgm = FT(1) - FT(1) / vgn
    θr = FT(0.075)
    ν = FT(0.287)
    θl_0 = FT(0.1)
    θl_surf = FT(0.267)
    Δt = FT(0.25)
    tf = FT(60 * 60 * 0.8)
    t0 = FT(0)

    msp = SoilParams{FT}(
        ν,
        vgn,
        vgα,
        vgm,
        ksat,
        θr,
        FT(1e-3),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )
    bottom_bc = FreeDrainage()
    top_bc = θDirichlet(θl_surf)
    domain = Domains.IntervalDomain(z₀, z₁, x3boundary = (:bottom, :top))
    mesh = Meshes.IntervalMesh(domain, nelems = n)

    cs = Spaces.CenterFiniteDifferenceSpace(mesh)
    fs = Spaces.FaceFiniteDifferenceSpace(cs)
    zc = Fields.coordinate_field(cs)


    θl = Fields.zeros(FT, cs) .+ θl_0

    # Solve Richard's Equation: ∂_t θl = ∂_z [K(θl)(∂_z ψ(θl) +1)]

    p = [zc, msp, top_bc, bottom_bc]
    function ∑tendencies!(dθl, θl, p, t)
        top = p[3]
        bot = p[4]
        compute_richards_rhs!(dθl, θl, p, top, bot)
    end
    #@show ∑tendencies!(similar(θl), θl, nothing, 0.0);

    # Solve the ODE operator

    prob = ODEProblem(∑tendencies!, θl, (t0, tf), p)
    sol = solve(
        prob,
        TsitPap8(),
        dt = Δt,
        saveat = 60 * Δt,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
    )

    bonan_sand_dataset = ArtifactWrapper(
        @__DIR__,
        isempty(get(ENV, "CI", "")),
        "richards_sand",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/2vk7bvyjah8xd5b7wxcqy72yfd2myjss.csv",
            filename = "sand_bonan_sp801.csv",
        ),],
    )
    datapath = get_data_folder(bonan_sand_dataset)
    data = joinpath(datapath, "sand_bonan_sp801.csv")
    ds_bonan = readdlm(data, ',')
    bonan_moisture = reverse(ds_bonan[:, 1])
    bonan_z = reverse(ds_bonan[:, 2]) ./ 100.0

    #    dirname = "richards"
    #    path = joinpath(@__DIR__, "test", "SoilModel", dirname)
    #    mkpath(path)
    #    plot(bonan_moisture, bonan_z, label = "Bonan solution", lw = 2, lc = :green)    
    #    plot!(
    #        parent(sol.u[end]),parent(zc),
    #        xlim = (0, 0.287),
    #        ylim = (-1.5, 0),
    #        label = "Clima solution",
    #        lc = :black,
    #        lw = 2,
    #        ls = :dash,
    #        legend = :outerright,
    #        xlabel = "θ(z)",
    #        ylabel = "z",
    #    )
    #    savefig(joinpath(path, "richards_sand.png"))
    @test sqrt.(sum((bonan_moisture .- parent(sol.u[end])) .^ 2.0)) < FT(0.1)
end
