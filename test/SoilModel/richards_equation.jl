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
    soil_rhs! = make_rhs(soil_model)

    land_model = LandHydrologyModel(soil = soil_model)
    Y2 = set_initial_state(land_model, (soil = initial_conditions,), 0.0)
    land_rhs! = make_rhs(land_model)
    dY1 = similar(Y1)
    dY2 = similar(Y2)
    a = soil_rhs!(dY1, Y1, [], 0.0)
    b = land_rhs!(dY2, Y2, [], 0.0)
    @test parent(Y1.soil.ϑ_l) ≈ parent(Y2.soil.ϑ_l)
    @test parent(b.soil.ϑ_l) ≈ parent(a.soil.ϑ_l)

end


@testset "Variably saturated equilibrium" begin
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
        θ_l = 0.494
        ρc_ds = model.soil_param_set.ρc_ds
        ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρc_ds, model.earth_param_set)
        ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, model.earth_param_set)
        return (ϑ_l = θ_l, θ_i = θ_i, ρe_int = ρe_int)
    end
    Y = set_initial_state(soil_model, initial_conditions, 0.0)
    soil_rhs! = make_rhs(soil_model)
    prob = ODEProblem(soil_rhs!, Y, (t0, tf), [])

    # solve simulation
    sol = solve(
        prob,
        SSPRK33(),
        dt = dt,
        saveat = 60 * dt,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
    )

    space_c, _ = make_function_space(domain)
    zc = Fields.coordinate_field(space_c)
    z = parent(zc)
    ϑ_l = [parent(sol.u[k].soil.ϑ_l) for k in 1:length(sol.u)]
    function expected(z, z_interface)
        ν = 0.495
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



    @test sqrt(mean(ϑ_l[end] .- expected.(z, -0.56)) .^ 2.0) < 1e-4

    #= 
    #Plotting
    indices = [1, 88, length(sol.t)]
    labels = ["IC", "6d", "36d"]
    plot1 = plot(
        xlim = (0.4, 0.525),
        ylim = (-10, 0),
        legend = :outerright,
        xlabel = "θ(z)",
        ylabel = "z",
    )
    for i in 1:1:length(indices)
        plot!(ϑ_l[indices[i]], z, label = labels[i], lw = 2)
    end

    plot!(expected.(z, -0.56), z, lw =1, label = "expected");
    plot(plot1)
        =#
end


@testset "Richards equation in sand, alternate BC" begin
    FT = Float64

    # General soil composition
    ν = FT(0.287)
    #Water specific
    Ksat = FT(34 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(3.96)
    vg_α = FT(2.7) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0.075)

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
    tf = FT(60 * 60 * 0.8)
    dt = FT(0.25)
    n = 150

    zmax = FT(0)
    zmin = FT(-1.5)
    domain = Column(FT, zlim = (zmin, zmax), nelements = n)

    #Boundary conditions
    top_state = (t) -> eltype(t)(0.267)
    bc = SoilDomainBC(
        domain;
        top = SoilComponentBC(hydrology = Dirichlet(top_state)),
        bottom = SoilComponentBC(hydrology = FreeDrainage()),
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

    function ic(z, t0, model)
        T = model.energy_model.T_profile(z, t0) # to be consistent with PrescribedT Default. 
        θ_i = 0.0
        θ_l = 0.1
        ρc_ds = model.soil_param_set.ρc_ds
        ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρc_ds, model.earth_param_set)
        ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, model.earth_param_set)
        return (ϑ_l = θ_l, θ_i = θ_i, ρe_int = ρe_int)
    end

    Y = set_initial_state(soil_model, ic, t0)
    soil_rhs! = make_rhs(soil_model)
    prob = ODEProblem(soil_rhs!, Y, (t0, tf), [])

    # solve simulation
    sol = solve(
        prob,
        SSPRK33(),
        dt = dt,
        saveat = 60 * dt,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
    )

    space_c, _ = make_function_space(domain)
    zc = Fields.coordinate_field(space_c)
    z = parent(zc)
    ϑ_l = parent(sol.u[193].soil.ϑ_l)


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
    @test sqrt.(sum((bonan_moisture .- ϑ_l) .^ 2.0)) < FT(0.1)
end
