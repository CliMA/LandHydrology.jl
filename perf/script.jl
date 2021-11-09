using ClimaCore: Fields
using CLIMAParameters
using CLIMAParameters.Planet: T_0
using UnPack
using OrdinaryDiffEq:
    ODEProblem,
    solve,
    CarpenterKennedy2N54,# does not work
    SSPRK33,
    SSPRK73
using LandHydrology
using LandHydrology.Models: default_initial_conditions
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
using LandHydrology.Domains: Column, make_function_space
using LandHydrology.SoilInterface
using LandHydrology.SoilInterface.SoilWaterParameterizations
using LandHydrology.SoilInterface.SoilHeatParameterizations
using LandHydrology.Simulations
using Statistics
using Test
using ArtifactWrappers
using DelimitedFiles



import Coverage
import Profile
const FT = Float64
function equilibrium_integrator()
    ν = FT(0.5)
    Ksat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    θ_r = FT(0)
    ν_ss_quartz = FT(0.92)
    ν_ss_om = FT(0.0)
    ν_ss_gravel = FT(0.0)
    κ_quartz = FT(7.7) # W/m/K
    κ_minerals = FT(2.5) # W/m/K
    κ_om = FT(0.25) # W/m/K
    κ_liq = FT(0.57) # W/m/K
    κ_ice = FT(2.29) # W/m/K
    κ_solid = k_solid(ν_ss_om, ν_ss_quartz, κ_quartz, κ_minerals, κ_om)
    κ_sat_frozen = ksat_frozen(κ_solid, ν, κ_ice)
    κ_sat_unfrozen = ksat_unfrozen(κ_solid, ν, κ_liq)
    ρc_ds = FT((1 - ν) * 1.926e06)
    #collect all params
    msp = SoilParams{FT}(
        ν = ν,
        S_s = S_s,
        ν_ss_gravel = ν_ss_gravel,
        ν_ss_om = ν_ss_om,
        ν_ss_quartz = ν_ss_quartz,
        ρc_ds = ρc_ds,
        κ_solid = κ_solid,
        κ_sat_unfrozen = κ_sat_unfrozen,
        κ_sat_frozen = κ_sat_frozen,
    )


    #Simulation and domain info
    t0 = FT(0)
    tf = FT(60 * 60 * 24 * 32)
    dt = FT(20)
    n = 20

    zmax = FT(0)
    zmin = FT(-2)
    domain = Column(FT, zlim = (zmin, zmax), nelements = n)

    #Boundary conditions
    top_flux = FT(0)
    bottom_flux = FT(0)
    bc = SoilColumnBC(;
        top = SoilComponentBC(
            hydrology = VerticalFlux(top_flux),
            energy = VerticalFlux(top_flux),
        ),
        bottom = SoilComponentBC(
            hydrology = VerticalFlux(bottom_flux),
            energy = VerticalFlux(bottom_flux),
        ),
    )

    # create model
    hydraulics_model =
        vanGenuchten{FT}(n = vg_n, α = vg_α, Ksat = Ksat, θr = θ_r)

    soil_model = SoilModel(
        FT;
        domain = domain,
        energy_model = SoilEnergyModel(),
        hydrology_model = SoilHydrologyModel{FT}(
            hydraulic_model = hydraulics_model,
        ),
        boundary_conditions = bc,
        soil_param_set = msp,
        earth_param_set = param_set,
    )

    # initial conditions
    function initial_conditions(z::FT, model::SoilModel)
        param_set = model.earth_param_set
        T = 289.0 + 5.0 * z
        θ_i = 0.0
        θ_l = 0.495
        ρcds = model.soil_param_set.ρc_ds
        ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρcds, param_set)
        ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, param_set)
        return (ϑ_l = θ_l, θ_i = θ_i, ρe_int = ρe_int)
    end
    Y, Ya = initialize_states(soil_model, initial_conditions, t0)
    soil_sim = Simulation(
        soil_model,
        SSPRK33(),
        Y_init = Y,
        dt = dt,
        tspan = (t0, tf),
        Ya_init = Ya,
    )
    return soil_sim.integrator
end

integrator = equilibrium_integrator()

step!(integrator, 1) # compile first
Profile.clear_malloc_data()
step!(integrator, 1)
#allocs = Coverage.analyze_malloc(".")
#=
julia> @benchmark step!(integrator,10.0)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  256.661 μs … 19.955 ms  ┊ GC (min … max):  0.00% … 97.76%
 Time  (median):     294.525 μs              ┊ GC (median):     0.00%
 Time  (mean ± σ):   381.341 μs ±  1.249 ms  ┊ GC (mean ± σ):  24.21% ±  7.23%

   ▄▁        █▆                                                 
  ▄████▄▃▃▃▂▅███▇▄▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂▂▂▂▂▂▂▁▂▂ ▃
  257 μs          Histogram: frequency by time          459 μs <

 Memory estimate: 758.80 KiB, allocs estimate: 4608.
=#
