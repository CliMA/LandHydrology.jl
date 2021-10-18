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
    hydraulics_model =
        vanGenuchten{FT}(n = vg_n, α = vg_α, Ksat = Ksat, θr = θ_r)

    soil_model = SoilModel(
        domain = domain,
        energy_model = PrescribedTemperatureModel(),
        hydrology_model = SoilHydrologyModel{FT}(
            hydraulic_model = hydraulics_model,
        ),
        boundary_conditions = bc,
        soil_param_set = msp,
        earth_param_set = param_set,
    )

    # initial conditions
    function initial_conditions(z, model)
        θ_i = 0.0
        θ_l = 0.494
        return (ϑ_l = θ_l, θ_i = θ_i)
    end
    Y, Ya = initialize_states(soil_model, initial_conditions, t0)
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

    z = parent(Ya.zc)
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
    hydraulics_model =
        vanGenuchten{FT}(n = vg_n, α = vg_α, Ksat = Ksat, θr = θ_r)
    # create model
    soil_model = SoilModel(
        domain = domain,
        energy_model = PrescribedTemperatureModel(),
        hydrology_model = SoilHydrologyModel{FT}(
            hydraulic_model = hydraulics_model,
        ),
        boundary_conditions = bc,
        soil_param_set = msp,
        earth_param_set = param_set,
    )

    # initial conditions

    function ic(z, model)
        θ_i = 0.0
        θ_l = 0.1
        return (ϑ_l = θ_l, θ_i = θ_i)
    end
    Y, Ya = initialize_states(soil_model, ic, t0)
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


    z = parent(Ya.zc)
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
