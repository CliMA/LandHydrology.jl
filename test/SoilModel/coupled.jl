@testset "Variably saturated equilibrium" begin
    FT = Float64
    ν = FT(0.5)
    Ksat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    θ_r = FT(0)
    ν_ss_quartz = FT(0.92)
    ν_ss_om = FT(0.0)
    ν_ss_gravel = FT(0.0)
    ρp = FT(2700)
    κ_quartz = FT(7.7) # W/m/K
    κ_minerals = FT(2.5) # W/m/K
    κ_om = FT(0.25) # W/m/K
    κ_liq = FT(0.57) # W/m/K
    κ_ice = FT(2.29) # W/m/K
    κ_solid = k_solid(ν_ss_om, ν_ss_quartz, κ_quartz, κ_minerals, κ_om)
    κ_sat_frozen = ksat_frozen(κ_solid, ν, κ_ice)
    κ_sat_unfrozen = ksat_unfrozen(κ_solid, ν, κ_liq)
    ρc_ds = FT((1 - ν) * 1.926e06)
    a = FT(0.24)
    b = FT(18.1)
    κ_dry_parameter = FT(0.053)
    #collect all params
    msp = SoilParams{FT}(
        ν,
        S_s,
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
    bc = SoilDomainBC(
        domain;
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
    soil_rhs! = make_rhs(soil_model)
    prob = ODEProblem(soil_rhs!, Y, (t0, tf), Ya)

    # solve simulation
    sol = solve(prob, SSPRK33(), dt = dt, saveat = 60 * dt)

    z = parent(Ya.zc)
    vlf = parent(sol.u[end].soil.ϑ_l)
    ρeint = parent(sol.u[end].soil.ρe_int)
    ρc_s = volumetric_heat_capacity.(vlf, 0.0, ρc_ds, Ref(param_set))

    temp = temperature_from_ρe_int.(ρeint, 0.0, ρc_s, Ref(param_set))
    function expected(z, z_interface)
        ν = 0.5
        S_s = 1e-3
        α = 2.6
        n = 2.0
        m = 0.5
        if z < z_interface
            return -S_s * (z - z_interface) + ν
        else
            return ν * (1 + (α * (z - z_interface))^n)^(-m)
        end
    end



    @test sqrt(mean(vlf .- expected.(z, -0.3)) .^ 2.0) < 1e-3
    @test sqrt(mean(temp .- 284.0) .^ 2.0) < 1e-3

end


@testset "test default ic" begin
    FT = Float64
    ν = FT(0.5)
    Ksat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    θ_r = FT(0)
    ν_ss_quartz = FT(0.92)
    ν_ss_om = FT(0.0)
    ν_ss_gravel = FT(0.0)
    ρp = FT(2700)
    κ_quartz = FT(7.7) # W/m/K
    κ_minerals = FT(2.5) # W/m/K
    κ_om = FT(0.25) # W/m/K
    κ_liq = FT(0.57) # W/m/K
    κ_ice = FT(2.29) # W/m/K
    κ_solid = k_solid(ν_ss_om, ν_ss_quartz, κ_quartz, κ_minerals, κ_om)
    κ_sat_frozen = ksat_frozen(κ_solid, ν, κ_ice)
    κ_sat_unfrozen = ksat_unfrozen(κ_solid, ν, κ_liq)
    ρc_ds = FT((1 - ν) * 1.926e06)
    a = FT(0.24)
    b = FT(18.1)
    κ_dry_parameter = FT(0.053)
    #collect all params
    msp = SoilParams{FT}(
        ν,
        S_s,
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
    bc = SoilDomainBC(
        domain;
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
        domain = domain,
        energy_model = SoilEnergyModel(),
        hydrology_model = SoilHydrologyModel{FT}(
            hydraulic_model = hydraulics_model,
        ),
        boundary_conditions = bc,
        soil_param_set = msp,
        earth_param_set = param_set,
    )

    Y_init = default_initial_conditions(soil_model)
    @test parent(Y_init.ϑ_l)[:] ≈ zeros(20) .+ 0.25
    @test parent(Y_init.θ_i)[:] ≈ zeros(20) .+ 0.0
    T0 = FT(T_0(soil_model.earth_param_set))

    ρc_s =
        volumetric_heat_capacity.(
            Y_init.ϑ_l,
            Y_init.θ_i,
            soil_model.soil_param_set.ρc_ds,
            Ref(soil_model.earth_param_set),
        )
    ρe_int =
        volumetric_internal_energy.(
            Y_init.θ_i,
            ρc_s,
            T0,
            Ref(soil_model.earth_param_set),
        )
    @test parent(Y_init.ρe_int)[:] ≈ parent(ρe_int)[:]


end
