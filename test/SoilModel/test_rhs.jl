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
        FT;
        domain = domain,
        energy_model = PrescribedTemperatureModel(T_profile = Tp),
        hydrology_model = PrescribedHydrologyModel(
            ϑ_l_profile = ϑ_lp,
            θ_i_profile = θ_ip,
        ),
        boundary_conditions = nothing,
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


# add other rhs tests here
