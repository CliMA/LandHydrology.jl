@testset "Turbulent surface fluxes" begin
    using SurfaceFluxes: DGScheme, surface_conditions
    using CLIMAParameters.Planet:
        R_v, grav, ρ_cloud_liq, R_d, cp_d, T_0, LH_v0, cp_v
    using Thermodynamics: Liquid, q_vap_saturation_generic, PhasePartition, cp_m

    FT = Float64


    ν = FT(0.55)
    vg_n = FT(1.68)
    vg_α = FT(5.0) # inverse meters
    θ_r = FT(0.084)

    #collect all params. We still dont have a good default for this.
    msp = SoilParams{FT}(ν = ν, ρc_ds = 1.0)

    n = 10

    zmax = FT(0)
    zmin = FT(-0.55)
    domain = Column(FT, zlim = (zmin, zmax), nelements = n)

    #Boundary conditions

    hm = vanGenuchten{FT}(; n = vg_n, α = vg_α, Ksat = 0.0, θr = θ_r)
    T_surf = 299.0
    ρ_a_sfc = 1.17
    q_atm = q_vap_saturation_generic(param_set, T_surf, ρ_a_sfc, Liquid())

    z_in = 0.05
    u_atm = 0.34
    surface_bc = PrescribedAtmosForcing{FT}(
        u_atm = u_atm,
        θ_atm = T_surf,
        z_atm = z_in,
        θ_scale = T_surf,
        ρ_a_sfc = ρ_a_sfc,
        q_atm = q_atm,
    )
    bc = SoilColumnBC(;
        top = surface_bc,
        bottom = SoilComponentBC(
            energy = VerticalFlux(0.0),
            hydrology = VerticalFlux(0.0),
        ),
    )
    soil_model = SoilModel(
        FT;
        domain = domain,
        energy_model = SoilEnergyModel(),
        hydrology_model = SoilHydrologyModel{FT}(hydraulic_model = hm),
        boundary_conditions = bc,
        soil_param_set = msp,
        earth_param_set = param_set,
    )

    # Test actual function used in RHS BC creation gives zero flux in case where we expect it to
    # initial conditions
    function initial_conditions(
        z::ft,
        model::SoilModel,
    ) where {ft <: AbstractFloat}
        param_set = model.earth_param_set
        ρcds = model.soil_param_set.ρc_ds
        ρc_s = volumetric_heat_capacity(
            model.soil_param_set.ν,
            ft(0.0),
            ρcds,
            param_set,
        )
        ρe_int = volumetric_internal_energy(ft(0.0), ρc_s, ft(299.0), param_set)
        return (ϑ_l = model.soil_param_set.ν, θ_i = ft(0.0), ρe_int = ρe_int)
    end
    Y, Ya = initialize_states(soil_model, initial_conditions, 0.0)
    soil_rhs! = make_rhs(soil_model)
    dY = similar(Y)
    soil_rhs!(dY, Y, Ya, 0.0)
    @test sum(parent(dY)) == 0.0

    ϑ_l = FT.([ν, ν + 1e-3, ν - 1e-3, ν])
    θ_i = FT.([0.0, 0.0, 0.0, 0.1])
    T = FT.([T_surf, T_surf, 289.5, 289.5])
    fluxes =
        compute_turbulent_surface_fluxes.(
            Ref(soil_model.energy_model),
            Ref(soil_model.hydrology_model),
            Ref(soil_model),
            ϑ_l,
            θ_i,
            T,
        )
    q_sat = q_vap_saturation_generic.(Ref(param_set), T, ρ_a_sfc, Ref(Liquid()))
    Lmo_guess = FT.([100.0 * z_in, u_atm, T_surf, q_atm])

    g = grav(param_set)
    Rv = R_v(param_set)

    correction =
        FT.([
            1.0,
            1.0,
            exp(
                g * matric_potential(hm, (ϑ_l[3] - hm.θr) / (ν - hm.θr)) / Rv /
                T[3],
            ),
            exp(
                g * matric_potential(
                    hm,
                    (ν - θ_i[4] - hm.θr) / (ν - θ_i[4] - hm.θr),
                ) / Rv / T[4],
            ),
        ])
    q_surf = correction .* q_sat

    _T_ref = FT(T_0(param_set))
    h_d = cp_d(param_set) .* (T .- _T_ref) .+ R_d(param_set) .* _T_ref
    lh = FT.(cp_v(param_set) .* (T .- _T_ref) .+ LH_v0(param_set))


    # Compute surface fluxes here and compare to result from function used within turbulent BC
    for i in [1, 2, 3, 4]
        conditions = surface_conditions(
            param_set,
            Lmo_guess,
            [u_atm, T_surf, q_atm],
            FT.([0.0, T[i], q_surf[i]]),
            FT.([0.001, 0.001, 0.001]),
            T_surf,
            z_in,
            DGScheme(),
        )
        ustar = conditions.x_star[1]
        tstar = conditions.x_star[2]
        qstar = conditions.x_star[3]
        qp = PhasePartition(q_surf[i])
        cpm = cp_m(param_set, qp)

        E = -ρ_a_sfc * ustar * qstar
        shf = -cpm * ρ_a_sfc * ustar * tstar .- h_d[i] * E
        lhf = lh[i] * E
        Ẽ = E / ρ_cloud_liq(param_set)

        heat_flux = shf + lhf
        @test fluxes[i][1] == heat_flux
        @test fluxes[i][2] == Ẽ
        if i == 2
            # Test that there is no shf if temps are the same
            @test tstar == 0.0
        end

    end

    # Test that oversaturated case gives same result as exactly saturated case
    @test fluxes[1] == fluxes[2]

    #test typing
    @test typeof(fluxes[1][1]) == FT
    @test typeof(fluxes[1][2]) == FT

    @test_throws Exception compute_turbulent_surface_fluxes(
        PrescribedTemperatureModel(),
        PrescribedHydrologyModel(),
        soil_model,
        ϑ_l[1],
        θ_i[1],
        T[1],
    )
    @test_throws Exception compute_turbulent_surface_fluxes(
        SoilEnergyModel(),
        PrescribedHydrologyModel(),
        soil_model,
        ϑ_l[1],
        θ_i[1],
        T[1],
    )
    @test_throws Exception compute_turbulent_surface_fluxes(
        PrescribedTemperatureModel(),
        SoilHydrologyModel{FT}(),
        soil_model,
        ϑ_l[1],
        θ_i[1],
        T[1],
    )


    @test_throws Exception boundary_fluxes(
        nothing,
        bc.top,
        :bottom,
        soil_model,
        nothing,
        nothing,
    )

end
