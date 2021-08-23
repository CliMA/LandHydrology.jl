import ClimaCore:
    Fields,
    Domains,
    Topologies,
    Meshes,
    DataLayouts,
    Operators,
    Geometry,
    Spaces
#using DiffEqCallbacks
using OrdinaryDiffEq: ODEProblem, solve, CarpenterKennedy2N54
using DelimitedFiles
using UnPack
using LandHydrology
using LandHydrology.SoilInterface
using LandHydrology.SoilWaterParameterizations
using ArtifactWrappers
using Test

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
const FT = Float64


abstract type bc end

struct θDirichlet{FT <: AbstractFloat} <: bc
    θvalue::FT
end
struct FreeDrainage <: bc end



function compute_richards_rhs!(dY, Y, p, top::θDirichlet,bottom::FreeDrainage)
    # θ_top = something, K∇h_bottom = K_bottom (∇h = 1). If ∇h = 1, θ_b = θ_f, so K = Kbottom at the face.
    (Y_hyd,) = Y.x
    (dY_hyd,) = dY.x
    @unpack    ϑ_l, θ_i = Y_hyd
    dϑ_l = dY_hyd.ϑ_l
    dθ_i = dY_hyd.θ_i
    sp = p[2]
    zc = p[1]
    @unpack ν, vgn, vgα, vgm, ksat, θr = sp
    θl = ϑ_l
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

    bc_b = Operators.SetValue(FT(1) * parent(K)[1])
    gradf2c = Operators.GradientF2C(bottom = bc_b) # set value on K∇h at bottom
    @. dϑ_l = gradf2c( If(K) * gradc2f(h))
    cs = axes(θ_i)
    dθ_i = Fields.zeros(eltype(θ_i),cs)
    return dY
    
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

    msp = SoilWaterParams{FT}(ν, vgn, vgα, vgm, ksat, θr, FT(1e-3))
    bottom_bc = FreeDrainage()
    top_bc = θDirichlet(θl_surf)
    domain = Domains.IntervalDomain(z₀, z₁, x3boundary = (:bottom, :top))
    mesh = Meshes.IntervalMesh(domain, nelems = n)

    cs = Spaces.CenterFiniteDifferenceSpace(mesh)
    fs = Spaces.FaceFiniteDifferenceSpace(cs)
    zc = Fields.coordinate_field(cs)
    
    function init_centers(zc)
        initial_value = 0.1
        θ_i = 0.0
        θl =  initial_value
        return (ϑ_l = θl, θ_i = θ_i,)
    end
    hydrology_model = SoilHydrologyModel(init_centers, nothing)
    energy_model = PrescribedTemperatureModel(nothing)
    soil_model = SoilModel(energy_model, hydrology_model, msp, param_set)
    Y = init_prognostic_vars(soil_model, cs)
    
    p = [zc, msp,top_bc, bottom_bc]
    function ∑tendencies!(dY, Y, p, t)
        top = p[3]
        bot = p[4]
        compute_richards_rhs!(dY, Y,p , top,bot)
    end
    
    prob = ODEProblem(∑tendencies!, Y, (t0, tf),p)
    sol = solve(
        prob,
        CarpenterKennedy2N54(),
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
    @test sqrt.(sum((bonan_moisture .- parent(sol.u[end].x[1].ϑ_l)).^2.0)) < FT(0.1)

end
