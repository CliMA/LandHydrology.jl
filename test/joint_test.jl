using LandHydrology.SurfaceFlowInterface

@testset "Nested soil rhs" begin
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

    ρc_ds = FT(2e6)
    #collect all params
    msp = SoilParams{FT}(
        ν,
        vg_n,
        vg_α,
        vg_m,
        Ksat,
        θ_r,
        S_s,
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
    bc = SoilDomainBC(
        domain;
        top = SoilComponentBC(hydrology = VerticalFlux(top_water_flux)),
        bottom = SoilComponentBC(hydrology = VerticalFlux(bottom_water_flux)),
    )

    # create model
    soil_model = SoilModel(
        domain = domain,
        energy_model = PrescribedTemperatureModel(),
        hydrology_model = SoilHydrologyModel(),
        boundary_conditions = bc,
        soil_param_set = msp,
        earth_param_set = param_set,
    )

    # initial conditions
    function initial_conditions(z, t0, model)
        T = model.energy_model.T_profile(z, t0) # to be consistent with PrescribedT Default. 
        θ_i = 0.0
        θ_l = 0.2 + 0.1 * z / 10.0
        ρc_ds = model.soil_param_set.ρc_ds
        ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρc_ds, model.earth_param_set)
        ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, model.earth_param_set)
        return (ϑ_l = θ_l, θ_i = θ_i, ρe_int = ρe_int)
    end
    Y1 = set_initial_state(soil_model, initial_conditions, 0.0)
    surface_model = SurfaceFlowModel()
    ic_surf = ()-> (h = [0.0],)
    land_model = LandHydrologyModel(soil = soil_model, surface_flow = surface_model)

    soil_rhs! = make_rhs(soil_model, land_model)
    Y2 = set_initial_state(land_model, (soil = initial_conditions,surface_flow = ic_surf), 0.0)
    land_rhs! = make_rhs(land_model)
    dY1 = similar(Y1)
    dY2 = similar(Y2)
    a = soil_rhs!(dY1, Y1, [], 0.0)
    b = land_rhs!(dY2, Y2, [], 0.0)
    @test parent(Y1.soil.ϑ_l) ≈ parent(Y2.soil.ϑ_l)
    @test parent(b.soil.ϑ_l) ≈ parent(a.soil.ϑ_l)


    t0 = FT(0)
    tf = FT(20)
    dt = FT(1)
    prob = ODEProblem(land_rhs!, Y2, (t0, tf), [])

    # solve simulation
    sol = solve(
        prob,
        SSPRK33(),
        dt = dt,
    )
    @test (sol.u[end].surface_flow.h[1] == 200.0)
end
