using ClimaCore: Fields
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using OrdinaryDiffEq:
    ODEProblem,
    solve,
    CarpenterKennedy2N54,# does not work
    SSPRK33
using LandHydrology
using LandHydrology.Domains
using LandHydrology.SoilInterface
using LandHydrology.SoilInterface.SoilHeatParameterizations
using Test
@testset "Prescribed 2x" begin
    FT = Float64
    ρc_ds = FT(2e6)
    msp = SoilParams{FT}(
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
        ρc_ds,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )
    t0 = FT(0)
    tf = FT(60)
    dt = FT(1)
    n = 50
    zmax = FT(0)
    zmin = FT(-10)
    domain = Column(FT, zlim = (zmin, zmax), nelements = n)
    soil_model = SoilModel(
        domain = domain,
        energy_model = PrescribedTemperatureModel(),
        hydrology_model = PrescribedHydrologyModel(ϑ_l_profile = (z,t)->eltype(z)(0.2), θ_i_profile = (z,t) -> eltype(z)(0.002*t)),
        soil_param_set = msp,
        earth_param_set = param_set,
    )

    function initial_conditions(z::Real, t0::Real, model::SoilModel)
        T = model.energy_model.T_profile(z, t0) # to be consistent with PrescribedT Default. 
        θ_i = 0.0
        θ_l = 0.2
        ρc_ds = model.soil_param_set.ρc_ds
        ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρc_ds, model.earth_param_set)
        ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, model.earth_param_set)
        return (ϑ_l = θ_l, θ_i = θ_i, ρe_int = ρe_int)
    end
    Y = set_initial_state(soil_model, initial_conditions, 0.0)
    soil_rhs! = make_rhs(soil_model)
    prob = ODEProblem(soil_rhs!, Y, (t0, tf), [])
    sol = solve(
        prob,
        SSPRK33(),
        dt = dt,
    )
    ϑ_l = sol.u[end].ϑ_l
    θ_i = sol.u[end].θ_i
    ρe_int = sol.u[end].ρe_int
    T = 288.0
    t = sol.t[end]
    θ_i_e = 0.002*t
    ϑ_l_e = 0.2
    ρc_ds = msp.ρc_ds
    ρc_s = volumetric_heat_capacity(ϑ_l_e, θ_i_e, ρc_ds, param_set)
    ρe_int_e = volumetric_internal_energy(θ_i_e, ρc_s, T, param_set)
    @test ρe_int .- ρe_int_e ≈ 0.0
    @test θ_i .- θ_i_e ≈ 0.0
    @test ϑ_l .- ϑ_l_e ≈ 0.0
end
