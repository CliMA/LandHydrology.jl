#  K = 0
# T < T_freeze in part of domain
# test that rhs = source
@testset "Phase Change" begin
    FT = Float64
    ν = FT(0.5)
    Ksat = FT(0)
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
    n = 20
    
    zmax = FT(0)
    zmin = FT(-2)
    Δz = FT((zmax-zmin)/n)
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
        sources = PhaseChange{FT}(Δz = Δz),
        soil_param_set = msp,
        earth_param_set = param_set,
    )
    
    # initial conditions
    function initial_conditions(z::FT, t0::FT, model::SoilModel)
        param_set = model.earth_param_set
        T = 270.0
        θ_i = 0.0
        θ_l = 0.4
        ρcds = model.soil_param_set.ρc_ds
        ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρcds, param_set)
        ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, param_set)
        return (ϑ_l = θ_l, θ_i = θ_i, ρe_int = ρe_int)
    end
    Y,Ya = set_initial_state(soil_model, initial_conditions, 0.0)
    soil_rhs! = make_rhs(soil_model)
    dY = similar(Y) .+ 0.0
    soil_rhs!(dY,Y,Ya,0.0)
    ρ_liq = ρ_cloud_liq(param_set)
    ρ_ice = ρ_cloud_ice(param_set)
    @test sum(parent(dY.soil.ϑ_l .*ρ_liq .+ dY.soil.θ_i .*ρ_ice) .== 0.0) .== 20

    computed_source! = add_source(PhaseChange{FT}(Δz = Δz), soil_model)
    dY2 = similar(Y) .+ 0.0
    computed_source!(dY2, Y, Ya, 0.0)
    @test parent(dY2.soil.ϑ_l) ≈ parent(dY.soil.ϑ_l)
    @test parent(dY2.soil.θ_i) ≈ parent(dY.soil.θ_i)
end

@testset "empty RHS unit test" begin
    FT = Float64
    n = 20
    zmax = FT(0)
    zmin = FT(-2)
    domain = Column(FT, zlim = (zmin, zmax), nelements = n)

    Tp(z, t) = 10.0 * z + t
    ϑ_lp(z, t) = 10.0 * z * t
    θ_ip(z, t) = 0.0

    soil_model = SoilModel(
        domain = domain,
        energy_model = PrescribedTemperatureModel(T_profile = Tp),
        hydrology_model = PrescribedHydrologyModel(
            ϑ_l_profile = ϑ_lp,
            θ_i_profile = θ_ip,
        ),
        boundary_conditions = nothing,
        soil_param_set = nothing,
        earth_param_set = nothing,
    )
    Ys = Dict()
    Y = Fields.FieldVector(; Ys...)
    t = 0.0
    space_c, _ = make_function_space(domain)
    zc = coordinates(space_c)
    p = initialize_auxiliary(soil_model, t, zc)
    soil_rhs! = make_rhs(soil_model)
    dY = similar(Y)
    soil_rhs!(dY, Y, p, t)
    @test dY == Y #dY still empty? not sure what a good test is.

    update_aux_en! = make_update_aux(soil_model.energy_model)
    update_aux_hydr! = make_update_aux(soil_model.hydrology_model)
    t = 10.0
    update_aux_en!(p, t)
    update_aux_hydr!(p, t)
    @test parent(p.soil.T) ≈ Tp.(parent(p.zc), t)
    @test parent(p.soil.ϑ_l) ≈ ϑ_lp.(parent(p.zc), t)
    @test parent(p.soil.θ_i) ≈ θ_ip.(parent(p.zc), t)

end
