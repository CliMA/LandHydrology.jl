#wow what a mess
import ClimaCore.Geometry
import ClimaCore:
    Fields,
    Domains,
    Topologies,
    Meshes,
    DataLayouts,
    Operators,
    Geometry,
    Spaces

using CLIMAParameters
using CLIMAParameters.Planet: ρ_cloud_liq, ρ_cloud_ice, cp_l, cp_i, T_0, LH_f0
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using RecursiveArrayTools
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33,Rosenbrock23, Tsit5,SSPRK432, Feagin14, TsitPap8,CarpenterKennedy2N54
using DelimitedFiles
using Statistics
using UnPack
using LandHydrology
using LandHydrology.SoilHeatParameterizations
using LandHydrology.SoilWaterParameterizations
using Test

const FT = Float64

@testset "Water and heat, dirichlet at top" begin
    ν = FT(0.395);
    ν_ss_quartz = FT(0.92)
    ν_ss_minerals = FT(0.08)
    ν_ss_om = FT(0.0)
    ν_ss_gravel = FT(0.0);
    
    Ksat = FT(4.42 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(1.89)
    vg_α = FT(7.5); # inverse meters
    vg_m = FT(1) -FT(1)/vg_n
    θ_r = FT(0)
    # Other constants needed:
    κ_quartz = FT(7.7) # W/m/K
    κ_minerals = FT(2.5) # W/m/K
    κ_om = FT(0.25) # W/m/K
    κ_liq = FT(0.57) # W/m/K
    κ_ice = FT(2.29); # W/m/K
    
    ρp = FT(2700); # kg/m^3
    κ_solid = k_solid(ν_ss_om, ν_ss_quartz, κ_quartz, κ_minerals, κ_om)
    κ_sat_frozen = ksat_frozen(κ_solid, ν, κ_ice)
    κ_sat_unfrozen = ksat_unfrozen(κ_solid, ν, κ_liq);
    
    ρc_ds = FT((1 - ν) * 1.926e06); # J/m^3/K
    a = FT(0.24)
    b = FT(18.1)
    κ_dry_parameter = FT(0.053)
    sp = SoilParams{FT}(ν,vg_n,vg_α,vg_m, Ksat, θ_r, S_s,
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
                        κ_dry_parameter)
    
    
    n = 50
    
    # Specify the domain boundaries
    zmax = FT(0)
    zmin = FT(-1)
    domain = Domains.IntervalDomain(zmin, zmax, x3boundary = (:bottom, :top))
    mesh = Meshes.IntervalMesh(domain, nelems = n)
    
    cs = Spaces.CenterFiniteDifferenceSpace(mesh)
    fs = Spaces.FaceFiniteDifferenceSpace(cs)
    zc = Fields.coordinate_field(cs)
    
    #initial conditions
    T_max = FT(289.0)
    T_min = FT(288.0)
    c = FT(20.0)
    T= @.  T_min + (T_max - T_min) * exp(-(zc - zmax) / (zmin - zmax) * c)
    
    θ_i = Fields.zeros(FT,cs)
    
    theta_max = FT(ν * 0.5)
    theta_min = FT(ν * 0.4)
    θ_l = @. theta_min + (theta_max - theta_min) * exp(-(zc - zmax) / (zmin - zmax) * c)
    
    
    ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
    ρe_int = volumetric_internal_energy.(θ_i, ρc_s, T, Ref(param_set))
    
    
    #BC
    Ttop = parent(T)[end]
    θ_top = parent(θ_l)[end]
    f_heat_bottom = FT(κ_sat_unfrozen)
    f_water_bottom = FT(Ksat)
    
    
    
    ##Prelim
    
    ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
    T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))
    κ_dry = k_dry(param_set, sp)
    S_r = relative_saturation.(θ_l, θ_i, ν)
    kersten = kersten_number.(θ_i, S_r, Ref(sp))
    κ_sat = saturated_thermal_conductivity.(
        θ_l,
        θ_i,
        κ_sat_unfrozen,
        κ_sat_frozen,
    )
    κ = thermal_conductivity.(κ_dry, kersten, κ_sat)
    
    S = effective_saturation.(θ_l; ν = ν, θr = θ_r)
    K = hydraulic_conductivity.(S; vgm = vg_m, ksat = Ksat)
    ψ = matric_potential.(S; vgn = vg_n, vgα = vg_α, vgm = vg_m)
    h = ψ .+ zc
    
    
    
    ρe_int_l = volumetric_internal_energy_liq.(T, Ref(param_set))
    
    
    ## Applying as dirichlet at the top, for water and temp, directly:
    
    gradc2f_heat = Operators.GradientC2F(top = Operators.SetValue(Ttop)) # for ∇T on face
    gradf2c_heat = Operators.GradientF2C(bottom = Operators.SetValue(f_heat_bottom))
    
    # Use state boundary values to get values of everything at the face - needed because the flux is not just ∇T, ∇h
    θ_i_top = parent(θ_i)[end]#we dont put BC on the ice...
    ρc_s_top = volumetric_heat_capacity(θ_top, θ_i_top, ρc_ds, param_set)
    S_r_top = relative_saturation(θ_top, θ_i_top, ν)
    kersten_top = kersten_number(θ_i_top, S_r_top, sp)
    κ_sat_top = saturated_thermal_conductivity(
        θ_top,
        θ_i_top,
        κ_sat_unfrozen,
        κ_sat_frozen,
    )
    κ_top = thermal_conductivity(κ_dry, kersten_top, κ_sat_top)
    
    
    
    ρe_int_l_top = volumetric_internal_energy_liq(Ttop, param_set)
    S_top = effective_saturation(θ_top; ν = ν, θr=θ_r)
    K_top = hydraulic_conductivity(S_top; vgm = vg_m, ksat = Ksat)
    ψ_top = matric_potential(S_top; vgn = vg_n, vgα = vg_α, vgm = vg_m)
    h_top = ψ_top
    
    #Interpolate κ and ρe_int_l *K to the face
    If_heat = Operators.InterpolateC2F(top = Operators.SetValue(κ_top))
    If_d_heat = Operators.InterpolateC2F(top = Operators.SetValue(K_top*ρe_int_l_top))
    
    
    gradc2f_water = Operators.GradientC2F(top = Operators.SetValue(h_top))# for ∇h on face
    
    dI_orig = @. gradf2c_heat(If_heat(κ)*gradc2f_heat(T) + If_d_heat(ρe_int_l*K)*gradc2f_water(h))
    

    ## now for water...
    If_water = Operators.InterpolateC2F(top = Operators.SetValue(K_top))
    dθ_orig  = @. Operators.GradientF2C(bottom = Operators.SetValue(f_water_bottom))( If_water(K)*gradc2f_water(h))
    
    

    ### vs estimating the fluxes
    # we need the top fluxes
    dz_top = cs.face_local_geometry.WJ[end]
    f_water_top = K_top*(h_top - parent(h)[end])/(dz_top*2)
    f_heat_top = κ_top*(Ttop - parent(T)[end])/(dz_top*2) .+ ρe_int_l_top*f_water_top
    
    interpc2f = Operators.InterpolateC2F()
    gradc2f_heat = Operators.GradientC2F()
    gradf2c_heat = Operators.GradientF2C(top = Operators.SetValue(f_heat_top), bottom = Operators.SetValue(f_heat_bottom))
    gradc2f_water = Operators.GradientC2F()
    gradf2c_water= Operators.GradientF2C(top = Operators.SetValue(f_water_top), bottom = Operators.SetValue(f_water_bottom))
    
    dI_flux = @. gradf2c_heat(interpc2f(κ) * gradc2f_heat(T) + interpc2f(ρe_int_l*K)*gradc2f_water(h))
    dθ_flux = @. gradf2c_water(interpc2f(K) * gradc2f_heat(h))
    
    @test sum(parent(dθ_orig) .== parent(dθ_flux)) .== n
    
    @test sum(parent(dI_orig) .== parent(dI_flux)) .== n
end


# Repeat but with dirichlet at top for heat and at bottom for water

@testset "Water and heat, dirichlet at different sides" begin
    #This is kind of a nightmare unless we want to make an approximation. see below ***
    ν = FT(0.395);
    ν_ss_quartz = FT(0.92)
    ν_ss_minerals = FT(0.08)
    ν_ss_om = FT(0.0)
    ν_ss_gravel = FT(0.0);
    
    Ksat = FT(4.42 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(1.89)
    vg_α = FT(7.5); # inverse meters
    vg_m = FT(1) -FT(1)/vg_n
    θ_r = FT(0)
    # Other constants needed:
    κ_quartz = FT(7.7) # W/m/K
    κ_minerals = FT(2.5) # W/m/K
    κ_om = FT(0.25) # W/m/K
    κ_liq = FT(0.57) # W/m/K
    κ_ice = FT(2.29); # W/m/K
    
    ρp = FT(2700); # kg/m^3
    κ_solid = k_solid(ν_ss_om, ν_ss_quartz, κ_quartz, κ_minerals, κ_om)
    κ_sat_frozen = ksat_frozen(κ_solid, ν, κ_ice)
    κ_sat_unfrozen = ksat_unfrozen(κ_solid, ν, κ_liq);
    
    ρc_ds = FT((1 - ν) * 1.926e06); # J/m^3/K
    a = FT(0.24)
    b = FT(18.1)
    κ_dry_parameter = FT(0.053)
    sp = SoilParams{FT}(ν,vg_n,vg_α,vg_m, Ksat, θ_r, S_s,
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
                        κ_dry_parameter)
    
    
    n = 50
    
    # Specify the domain boundaries
    zmax = FT(0)
    zmin = FT(-1)
    domain = Domains.IntervalDomain(zmin, zmax, x3boundary = (:bottom, :top))
    mesh = Meshes.IntervalMesh(domain, nelems = n)
    
    cs = Spaces.CenterFiniteDifferenceSpace(mesh)
    fs = Spaces.FaceFiniteDifferenceSpace(cs)
    zc = Fields.coordinate_field(cs)
    
    #initial conditions
    T_max = FT(289.0)
    T_min = FT(288.0)
    c = FT(20.0)
    T= @.  T_min + (T_max - T_min) * exp(-(zc - zmax) / (zmin - zmax) * c)
    
    θ_i = Fields.zeros(FT,cs)
    
    theta_max = FT(ν * 0.5)
    theta_min = FT(ν * 0.4)
    θ_l = @. theta_min + (theta_max - theta_min) * exp(-(zc - zmax) / (zmin - zmax) * c)
    
    
    ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
    ρe_int = volumetric_internal_energy.(θ_i, ρc_s, T, Ref(param_set))
    
    
    #BC
    Ttop = parent(T)[end]
    θ_bottom = parent(θ_l)[1]
    f_heat_bottom = FT(κ_sat_unfrozen)
    f_water_top = FT(Ksat)
    
    
    
    ##Prelim
   
    ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
    T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))
    κ_dry = k_dry(param_set, sp)
    S_r = relative_saturation.(θ_l, θ_i, ν)
    kersten = kersten_number.(θ_i, S_r, Ref(sp))
    κ_sat = saturated_thermal_conductivity.(
        θ_l,
        θ_i,
        κ_sat_unfrozen,
        κ_sat_frozen,
    )
    κ = thermal_conductivity.(κ_dry, kersten, κ_sat)
    
    S = effective_saturation.(θ_l; ν = ν, θr = θ_r)
    K = hydraulic_conductivity.(S; vgm = vg_m, ksat = Ksat)
    ψ = matric_potential.(S; vgn = vg_n, vgα = vg_α, vgm = vg_m)
    h = ψ .+ zc
    
    
    
    ρe_int_l = volumetric_internal_energy_liq.(T, Ref(param_set))
    
    
    ## Applying as dirichlet at the top, for Temp, and bottom, for water:
    
    gradc2f_heat = Operators.GradientC2F(top = Operators.SetValue(Ttop)) # for ∇T on face
    gradf2c_heat = Operators.GradientF2C(bottom = Operators.SetValue(f_heat_bottom))


    # Next up, we need κ at the face at the top. κ depends on water content, but at the top, we only know water flux.

    #*****issue: Kface*(hface-hcenter)/dz -> flux constraint at top implies a constraint on theta face! but it would be a pain to solve for it. nonlinear.
    θ_top = parent(θ_l)[end]## stand in to proceed.

    θ_i_top = parent(θ_i)[end]#we dont put BC on the ice...
    ρc_s_top = volumetric_heat_capacity(θ_top, θ_i_top, ρc_ds, param_set)
    S_r_top = relative_saturation(θ_top, θ_i_top, ν)
    kersten_top = kersten_number(θ_i_top, S_r_top, sp)
    κ_sat_top = saturated_thermal_conductivity(
        θ_top,
        θ_i_top,
        κ_sat_unfrozen,
        κ_sat_frozen,
    )
    κ_top = thermal_conductivity(κ_dry, kersten_top, κ_sat_top)
    
    
    
    ρe_int_l_top = volumetric_internal_energy_liq(Ttop, param_set)
    S_bottom = effective_saturation(θ_bottom; ν = ν, θr=θ_r)
    S_top = effective_saturation(θ_top; ν = ν, θr=θ_r)
    K_top = hydraulic_conductivity(S_top; vgm = vg_m, ksat = Ksat)
    #*****In full model, K depends on T as well - so to be fully consistent, we'd have to solve for T_face at the bottom resulting from the flux BC. Simpler thing is to use Extrapolate()
    K_bottom = hydraulic_conductivity(S_bottom; vgm = vg_m, ksat = Ksat)
    
    ψ_bottom = matric_potential(S_bottom; vgn = vg_n, vgα = vg_α, vgm = vg_m)
    h_bottom = ψ_bottom + zmin #bottom of domain
    
    #Interpolate κ and ρe_int_l *K to the face
    If_heat = Operators.InterpolateC2F(top = Operators.SetValue(κ_top))
    If_d_heat = Operators.InterpolateC2F(top = Operators.SetValue(K_top*ρe_int_l_top))
    
    
    gradc2f_water = Operators.GradientC2F(top= Operators.SetGradient(f_water_top/K_top))
    
    dI_orig = @. gradf2c_heat(If_heat(κ)*gradc2f_heat(T) + If_d_heat(ρe_int_l*K)*gradc2f_water(h)) # this gives correct water flux at top, and ρe_int is correct on top b/c we have T_top
    #*****only error ? is that κ is not consistent with water flux at the top.?
    

    ## now for water...
    gradc2f_water = Operators.GradientC2F(bottom= Operators.SetValue(h_bottom))# for ∇h on face
    
    If_water = Operators.InterpolateC2F(bottom = Operators.SetValue(K_bottom))
    dθ_orig  = @. Operators.GradientF2C(top = Operators.SetValue(f_water_top))( If_water(K)*gradc2f_water(h))
    
    

    ### vs estimating the fluxes
    # we need the top fluxes
    dz_top = cs.face_local_geometry.WJ[end]
    dz_bottom = cs.face_local_geometry.WJ[1]
    f_water_bottom = K_bottom*(parent(h)[1]-h_bottom)/(dz_bottom*2)
    f_heat_top = κ_top*(Ttop - parent(T)[end])/(dz_top*2) .+ ρe_int_l_top*f_water_top
    
    interpc2f = Operators.InterpolateC2F()
    gradc2f_heat = Operators.GradientC2F()
    gradf2c_heat = Operators.GradientF2C(top = Operators.SetValue(f_heat_top), bottom = Operators.SetValue(f_heat_bottom))
    gradc2f_water = Operators.GradientC2F()
    gradf2c_water= Operators.GradientF2C(top = Operators.SetValue(f_water_top), bottom = Operators.SetValue(f_water_bottom))
    
    dI_flux = @. gradf2c_heat(interpc2f(κ) * gradc2f_heat(T) + interpc2f(ρe_int_l*K)*gradc2f_water(h))
    dθ_flux = @. gradf2c_water(interpc2f(K) * gradc2f_heat(h))
    
    @test sum(parent(dθ_orig) .== parent(dθ_flux)) == n
    
    @test sum(parent(dI_orig) .== parent(dI_flux)) == n
end



@testset "Water alone, Dirichlet at top" begin
    n = 150
    z₀ = FT(-1.5)
    z₁ = FT(0)
    ksat = FT(34 / (3600 * 100))
    vgn = FT(3.96)
    vgα = FT(2.7)
    vgm = FT(1)- FT(1)/vgn
    θr = FT(0.075)
    ν = FT(0.287)
    θl_0 = FT(0.1)

    domain = Domains.IntervalDomain(z₀, z₁, x3boundary = (:bottom, :top))
    mesh = Meshes.IntervalMesh(domain, nelems = n)
    
    cs = Spaces.CenterFiniteDifferenceSpace(mesh)
    fs = Spaces.FaceFiniteDifferenceSpace(cs)
    zc = Fields.coordinate_field(cs)
    
    
    θl = Fields.zeros(FT, cs) .+ θl_0
    S = effective_saturation.(θl; ν = ν, θr = θr)
    K = hydraulic_conductivity.(S; vgm = vgm, ksat = ksat)
    ψ = matric_potential.(S; vgn = vgn, vgα = vgα, vgm = vgm)
    h = ψ .+ zc

    θ_top = FT(0.267)
    S_top = effective_saturation.(θ_top; ν = ν, θr = θr)
    h_top = matric_potential(S_top; vgn = vgn, vgα = vgα, vgm = vgm)#ztop = 0
    K_top = hydraulic_conductivity(S_top; vgm = vgm, ksat = ksat)
    bc_ktop = Operators.SetValue(K_top)
    #we'd do the following to get the space inside RHS
    center_space = axes(θl)
    dz_top = center_space.face_local_geometry.WJ[end]
    #check why factor of 2 shows up
    flux_top = K_top*(h_top - parent(h)[end])/(dz_top)/FT(2)
    flux_bottom = parent(K)[1]
    If = Operators.InterpolateC2F(;
                                  top = bc_ktop)
    

    #original
    gradc2f = Operators.GradientC2F(top  = Operators.SetValue(h_top))
    gradf2c = Operators.GradientF2C(bottom = Operators.SetValue(flux_bottom)) # set value on K∇h at bottom
    original = @. gradf2c( If(K) * gradc2f(h))

    #new
    gradc2f_new  = Operators.GradientC2F()
    If_new = Operators.InterpolateC2F()
    gradf2c_new = Operators.GradientF2C(top = Operators.SetValue(flux_top), bottom = Operators.SetValue(flux_bottom))
    new = @. gradf2c_new(If_new(K)*gradc2f_new(h))

    @test sum(parent(original) .== parent(new)) == n
    
end

@testset "Heat alone, Dirichlet at top/bottom" begin
    ν = 0.495
    ν_ss_gravel = 0.1
    ν_ss_om = 0.1
    ν_ss_quartz = 0.1
    ρc_ds = 0.43314518988433487
    κ_solid = 8.0
    ρp = 2700.0
    κ_sat_unfrozen = 0.57
    κ_sat_frozen = 2.29
    a = 0.24
    b = 18.1
    κ_dry_parameter = 0.053
    sp = SoilHeatParams(
        ν,
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
        κ_dry_parameter
    )
    n = 60
    
    # Specify the domain boundaries
    zmax = FT(1)
    zmin = FT(0)
    domain = Domains.IntervalDomain(zmin, zmax, x3boundary = (:bottom, :top))
    mesh = Meshes.IntervalMesh(domain, nelems = n)
    
    cs = Spaces.CenterFiniteDifferenceSpace(mesh)
    fs = Spaces.FaceFiniteDifferenceSpace(cs)
    zc = Fields.coordinate_field(cs)
    
    T = zc
    Ttop = FT(10.0)
    Tbottom = FT(-5.0)
    

    
    θ_l = Fields.zeros(FT,cs) .+0.2
    θ_i = Fields.zeros(FT,cs)
    κ_dry = k_dry(param_set, sp)
    S_r = relative_saturation.(θ_l, θ_i, ν)
    kersten = kersten_number.(θ_i, S_r, Ref(sp))
    κ_sat = saturated_thermal_conductivity.(
        θ_l,
        θ_i,
        κ_sat_unfrozen,
        κ_sat_frozen,
    )
    κ = thermal_conductivity.(κ_dry, kersten, κ_sat)
   
    gradc2f = Operators.GradientC2F(top = Operators.SetValue(Ttop), bottom = Operators.SetValue(Tbottom))
    gradf2c = Operators.GradientF2C()
    If = Operators.InterpolateC2F(top = Operators.Extrapolate(), bottom = Operators.Extrapolate())# not sure about this. look into case in integrated model.
    original =  @.  gradf2c(If(κ) * gradc2f(T))



    center_space = axes(T)
    dz_top = center_space.face_local_geometry.WJ[end]
    dz_bottom = center_space.face_local_geometry.WJ[1]
    #check why factor of 2 shows up
    flux_top = parent(κ)[end]*(Ttop - parent(T)[end])/(dz_top)/FT(2)
    flux_bottom = parent(κ)[1]*(parent(T)[1]- Tbottom)/(dz_bottom)/FT(2)




    
    gradc2f_new  = Operators.GradientC2F()
    If_new = Operators.InterpolateC2F()
    gradf2c_new = Operators.GradientF2C(top = Operators.SetValue(flux_top), bottom = Operators.SetValue(flux_bottom))
    new = @. gradf2c_new(If_new(κ)*gradc2f_new(T))

    @test sum(parent(original) .== parent(new)) == n

end
