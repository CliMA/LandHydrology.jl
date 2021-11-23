using ClimaCore: Fields
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
param_set = EarthParameterSet()
using CLIMAParameters.Planet: T_0
using UnPack
using OrdinaryDiffEq:
    ODEProblem,
    solve,
    SSPRK33,
    SSPRK73
using LandHydrology
using LandHydrology: make_rhs, LandHydrologyModel, initialize_land_states, default_land_initial_conditions
using LandHydrology.SubComponentModels: default_initial_conditions, make_tendency_terms, make_update_aux, NotIncluded
using LandHydrology.Domains: Column, make_function_space
using LandHydrology.SoilInterface
using LandHydrology.SoilInterface.SoilWaterParameterizations
using LandHydrology.SoilInterface.SoilHeatParameterizations
using LandHydrology.Simulations
using LandHydrology.SurfaceWater
using Statistics
using Test
using ArtifactWrappers
using DelimitedFiles
using LandHydrology: PrescribedAtmosState
using LandHydrology.SurfaceWater: SurfaceWaterModel

@testset "No surface water" begin
    FT = Float64
    
    # General soil composition
    ν = FT(0.495)
    #Water specific
    Ksat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0)
    
    #collect all params
    msp = SoilParams{FT}(ν = ν, S_s = S_s)
    
    
    #Simulation and domain info
    t0 = FT(0)
    tf = FT(60 * 60 * 24 * 36)
    dt = FT(100)
    n = 50
    
    zmax = FT(0)
    zmin = FT(-10)
    domain = Column(FT, zlim = (zmin, zmax), nelements = n)
    
    #Boundary conditions
    top_water_flux = FT(0)
    bottom_water_flux = FT(0)
    bc = SoilColumnBC(;
                      top = SoilComponentBC(hydrology = SurfaceWaterBC(),),
                      bottom = SoilComponentBC(hydrology = VerticalFlux(bottom_water_flux)),
                      )
    
    # create model
    hydraulics_model =
        vanGenuchten{FT}(n = vg_n, α = vg_α, Ksat = Ksat, θr = θ_r)
    
    soil_model = SoilModel(
        FT;
        domain = domain,
        energy_model = PrescribedTemperatureModel(),
        hydrology_model = SoilHydrologyModel{FT}(
            hydraulic_model = hydraulics_model,
        ),
        boundary_conditions = bc,
        soil_param_set = msp,
    )
    
    surface = NotIncluded()
    
    precip = FT(-5e-7)
    land_model = LandHydrologyModel{FT}(soil_model,surface, param_set; atmos_state = PrescribedAtmosState{FT}(precip))

    function initial_conditions(z, land_model, model)
        θ_i = 0.0
        θ_l = 0.4
        return (ϑ_l = θ_l, θ_i = θ_i)
    end
    Y, Ya = initialize_land_states(land_model, (;soil =initial_conditions))
    land_rhs! = make_rhs(land_model)
    land_sim = Simulation(
        land_model,
        SSPRK33(),
        Y_init = Y,
        dt = dt,
        tspan = (t0, tf),
        Ya_init = Ya,
        saveat = 60 * dt,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
    )
    
    # solve simulation
    @test step!(land_sim) isa Nothing # either error or integration runs
    run!(land_sim)
    sol = land_sim.integrator.sol
    
    z = parent(Ya.soil.zc)
    # seems reasonable
end


@testset "With surface water" begin
    FT = Float64
    
    # General soil composition
    ν = FT(0.495)
    #Water specific
    Ksat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0)
    
    #collect all params
    msp = SoilParams{FT}(ν = ν, S_s = S_s)
    
    
    #Simulation and domain info
    t0 = FT(0)
    tf = FT(60 * 60 * 24 * 36)
    dt = FT(100)
    n = 50
    
    zmax = FT(0)
    zmin = FT(-10)
    domain = Column(FT, zlim = (zmin, zmax), nelements = n)
    
    #Boundary conditions
    top_water_flux = FT(0)
    bottom_water_flux = FT(0)
    bc = SoilColumnBC(;
                      top = SoilComponentBC(hydrology = SurfaceWaterBC(),),
                      bottom = SoilComponentBC(hydrology = VerticalFlux(bottom_water_flux)),
                      )
    
    # create model
    hydraulics_model =
        vanGenuchten{FT}(n = vg_n, α = vg_α, Ksat = Ksat, θr = θ_r)
    
    soil_model = SoilModel(
        FT;
        domain = domain,
        energy_model = PrescribedTemperatureModel(),
        hydrology_model = SoilHydrologyModel{FT}(
            hydraulic_model = hydraulics_model,
        ),
        boundary_conditions = bc,
        soil_param_set = msp,
    )
    
    
    surface = SurfaceWaterModel{FT}()
    function initial_conditions_sfc()
        return (;h = 0.0)
    end
    
    precip = FT(-5e-7)
    land_model = LandHydrologyModel{FT}(soil_model,surface, param_set; atmos_state = PrescribedAtmosState{FT}(precip))
    function initial_conditions(z, _...)
        θ_i = 0.0
        θ_l = 0.4
        return (ϑ_l = θ_l, θ_i = θ_i)
    end
    Y, Ya = initialize_land_states(land_model, (;soil =initial_conditions,sfc_water = initial_conditions_sfc))
    land_rhs! = make_rhs(land_model)
    land_sim = Simulation(
        land_model,
        SSPRK33(),
        Y_init = Y,
        dt = dt,
        tspan = (t0, tf),
        Ya_init = Ya,
        saveat = 60 * dt,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
    )
    
    # solve simulation
    @test step!(land_sim) isa Nothing # either error or integration runs
    run!(land_sim)
    sol = land_sim.integrator.sol
    
    z = parent(Ya.soil.zc)
    h = [sol.u[k].sfc_water.h for k in 1:520]
    # seems reasonable
end

