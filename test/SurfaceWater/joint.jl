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


    # Define the soil model
    # Soil parameters - stored in the soil model structure
    ν = FT(0.495)
    Ksat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0)
    msp = SoilParams{FT}(ν = ν, S_s = S_s)
    
    #domain info for the soil model
    n = 50
    zmax = FT(0)
    zmin = FT(-10)
    domain = Column(FT, zlim = (zmin, zmax), nelements = n)
    
    # Boundary conditions for the soil - indicates that at the top of the domain,
    # we have a surface water boundary condition (ponding), while at the bottom,
    # we have an impermeable surface.
    bottom_water_flux = FT(0)
    bc = SoilColumnBC(;
                      top = SoilComponentBC(hydrology = SurfaceWaterBC(),),
                      bottom = SoilComponentBC(hydrology = VerticalFlux(bottom_water_flux)),
                      )
    
    # A parameterization choice for the soil hydraulics - the van Genuchten model
    # The user could also choose Brooks and Corey here.
    hydraulics_model =
        vanGenuchten{FT}(n = vg_n, α = vg_α, Ksat = Ksat, θr = θ_r)

    # Here is the soil model. We have to give a domain, information about the soil energy and hydrology
    # equations neededs, boundary conditions, and soil parameters.
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
    
    # The surface water is comparatively very simple. It is just a pond on top of the soil,
    # there are no parameters, or boundary conditions to give. In fact, for a column model,
    # we don't even need to supply a domain, but in the future we would (e.g. a location on the
    # surface where that surface water is).
    surface = SurfaceWaterModel{FT}()

    # In this case, we want to run the model with a prescribed atmos boundary condition (like reanalysis data)
    # Right now this only accepts `precip` as a field, but in the future this would also have wind speed, temperature,
    # specific humidity, etc as a function of time and space. This is both how we pass in the info to the RHS
    # and indicate what type of atmosphere the land will be interacting with. Other options would include NoAtmosphere
    # (for simple simulations of e.g. lab soil column tests) or a CoupledAtmosphere.
    precip = FT(-5e-7)
    atmos_state = PrescribedAtmosState{FT}(precip)
    # Now we make the land model
    # param_set has global earth parameters (e.g. density of water...)
    land_model = LandHydrologyModel{FT}(soil_model,
                                        surface,
                                        param_set;
                                        atmos_state =  atmos_state)
    # Next up, we need initial conditions for the LandModel.
    # These are only required for variables with a time derivative.
    # all other initial vaues are computed internally.
    # These are defined by model.
    
    function initial_conditions_sfc()
        return (;h = 0.0)
    end
    function initial_conditions(z, _...)
        θ_i = 0.0
        θ_l = 0.4
        return (ϑ_l = θ_l, θ_i = θ_i)
    end
    Y, Ya = initialize_land_states(land_model, (;soil =initial_conditions,sfc_water = initial_conditions_sfc))

    # Now we have the inital prognostic state Y and aux state Ya for the land model. 
    land_rhs! = make_rhs(land_model)
    # this function computes the entire right hand side of dY/dt = RHS(Y, Ya, t)

    # time step info
    t0 = FT(0)
    tf = FT(60 * 60 * 24 * 36)
    dt = FT(100)

    # simulation info
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

