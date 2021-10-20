@testset "Heat Analytic Unit Test" begin
    FT = Float64
    ν = FT(0.495)
    ν_ss_gravel = FT(0.1)
    ν_ss_om = FT(0.1)
    ν_ss_quartz = FT(0.1)
    ρc_ds = FT(0.43314518988433487)
    κ_solid = FT(8.0)
    ρp = FT(2700.0)
    κ_sat_unfrozen = FT(0.57)
    κ_sat_frozen = FT(2.29)
    a = FT(0.24)
    b = FT(18.1)
    κ_dry_parameter = FT(0.053)
    msp = SoilParams{FT}(
        ν,
        0.0,
        ν_ss_gravel,
        ν_ss_om,
        ν_ss_quartz,
        ρc_ds,
        κ_solid,
        ρp,
        κ_sat_unfrozen,
        κ_sat_frozen,
        a,
        b,
        κ_dry_parameter,
    )


    t0 = FT(0)
    tf = FT(2)
    dt = FT(1e-4)
    n = 60

    # Specify the domain boundaries
    zmax = FT(1)
    zmin = FT(0)
    domain = Column(FT, zlim = (zmin, zmax), nelements = n)

    tau = FT(1) # period (sec)
    A = FT(5) # amplitude (K)
    ω = FT(2 * pi / tau)
    topbc = Dirichlet(t -> typeof(t)(0.0))
    bottombc = Dirichlet(t -> typeof(t)(A * cos(ω * t)))
    bc = SoilDomainBC(
        domain;
        top = SoilComponentBC(energy = topbc),
        bottom = SoilComponentBC(energy = bottombc),
    )

    # create model
    soil_model = SoilModel(
        domain = domain,
        energy_model = SoilEnergyModel(),
        hydrology_model = PrescribedHydrologyModel(),
        boundary_conditions = bc,
        soil_param_set = msp,
        earth_param_set = param_set,
    )

    # initial conditions
    @test_throws ErrorException default_initial_conditions(soil_model)
    function energy_ic(z, model)
        T = 0.0
        θ_i = 0.0
        θ_l = 0.0
        ρc_s = volumetric_heat_capacity(
            θ_l,
            θ_i,
            model.soil_param_set.ρc_ds,
            model.earth_param_set,
        )
        ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, model.earth_param_set)
        return (; ρe_int = ρe_int)
    end
    Y, Ya = initialize_states(soil_model, energy_ic, t0)
    soil_sim = Simulation(
        soil_model,
        SSPRK33(),
        Y_init = Y,
        dt = dt,
        tspan = (t0, tf),
        p = Ya,
        saveat = 60 * dt,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
    )

    # solve simulation
    @test step!(soil_sim) isa Nothing # either error or integration runs
    run!(soil_sim)
    sol = soil_sim.integrator.sol
    t = sol.t
    z = parent(Ya.zc)[:]
    num =
        exp.(sqrt(ω / 2) * (1 + im) * (1 .- z)) .-
        exp.(-sqrt(ω / 2) * (1 + im) * (1 .- z))
    denom = exp(sqrt(ω / 2) * (1 + im)) - exp.(-sqrt(ω / 2) * (1 + im))
    analytic_soln = real(num .* A * exp(im * ω * tf) / denom)
    ρe_intf = parent(sol.u[end].soil.ρe_int)[:]
    θ_lf = soil_model.hydrology_model.ϑ_l_profile.(z, tf)
    θ_if = soil_model.hydrology_model.ϑ_l_profile.(z, tf)
    ρc_sf = volumetric_heat_capacity.(θ_lf, θ_if, ρc_ds, Ref(param_set))
    Tfinal = temperature_from_ρe_int.(ρe_intf, θ_if, ρc_sf, Ref(param_set))
    MSE = mean((analytic_soln .- Tfinal) .^ 2.0)
    @test MSE < 1e-6
end
