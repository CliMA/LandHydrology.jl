import ClimaCore.Geometry, LinearAlgebra, UnPack
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
using Plots
using DelimitedFiles
using UnPack
using LandHydrology
using LandHydrology.SoilHeatParameterizations
using LandHydrology.SoilWaterParameterizations
using Test
using Statistics

const FT = Float64


function compute_integrated_rhs!(dY, Y, t, p)

    sp = p[1]
    param_set = p[2]
    zc = p[3]
    @unpack ν,vgn,vgα,vgm,ksat,θr,ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp
    ϑ_l = Y.x[1]
    θ_i = Y.x[2]
    ρe_int = Y.x[3]
    f = eltype(θ_i)
    dϑ_l = dY.x[1]
    dθ_i = dY.x[2]
    dρe_int = dY.x[3]

    # Compute center values of everything
    θ_l = ϑ_l
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



    S = effective_saturation.(θ_l; ν = ν, θr = θr)
    K = hydraulic_conductivity.(S; vgm = vgm, ksat = ksat)
    ψ = matric_potential.(S; vgn = vgn, vgα = vgα, vgm = vgm)
    h = ψ .+ zc



    ρe_int_l = volumetric_internal_energy_liq.(T, Ref(param_set))
    # Boundary conditions on the Flux. being sloppy with negative sign.
    top_heat_flux = f(0.0) # on κ∇T
    top_water_flux = f(0.0)
    bottom_water_flux = f(0.0) # on K∇h
    bottom_heat_flux = f(0.0)

    interpc2f = Operators.InterpolateC2F()
    gradc2f_heat = Operators.GradientC2F()
    gradf2c_heat = Operators.GradientF2C(top = Operators.SetValue(top_heat_flux), bottom = Operators.SetValue(bottom_heat_flux))

    gradc2f_water = Operators.GradientC2F()
    gradf2c_water= Operators.GradientF2C(top = Operators.SetValue(top_water_flux), bottom = Operators.SetValue(bottom_water_flux))

    @. dϑ_l = gradf2c_water( interpc2f(K) * gradc2f_water(h)) #Richards equation
    @. dρe_int = gradf2c_heat(interpc2f(κ) * gradc2f_heat(T) + interpc2f(ρe_int_l*K)*gradc2f_water(h))
    @. dθ_i = K - K#/Fields.zeros(f,cs)

    return dY
  
end


ν = FT(0.395);
# Soil solids
# are the components of soil besides water, ice, gases, and air.
# We specify the soil component fractions, relative to all soil solids.
# These should sum to unity; they do not account for pore space.
ν_ss_quartz = FT(0.92)
ν_ss_minerals = FT(0.08)
ν_ss_om = FT(0.0)
ν_ss_gravel = FT(0.0);
# Other parameters include the hydraulic conductivity at saturation, the specific
# storage, and the van Genuchten parameters for sand.
# We recommend Chapter 8 of  [Bonan19a](@cite) for finding parameters
# for other soil types.
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
# The particle density of organic material-free soil is
# equal to the particle density of quartz and other minerals ([BallandArp2005](@cite)):
ρp = FT(2700); # kg/m^3
# We calculate the thermal conductivities for the solid material
# and for saturated soil. These functions are taken from [BallandArp2005](@cite).
κ_solid = k_solid(ν_ss_om, ν_ss_quartz, κ_quartz, κ_minerals, κ_om)
κ_sat_frozen = ksat_frozen(κ_solid, ν, κ_ice)
κ_sat_unfrozen = ksat_unfrozen(κ_solid, ν, κ_liq);
# Next, we calculate the volumetric heat capacity of dry soil. Dry soil
# refers to soil that has no water content.
ρc_ds = FT((1 - ν) * 1.926e06); # J/m^3/K
a = FT(0.24)
b = FT(18.1)
κ_dry_parameter = FT(0.053)
msp = SoilParams{FT}(ν,vg_n,vg_α,vg_m, Ksat, θ_r, S_s,
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



t0 = FT(0)
tf = FT(60 * 60 * 72)
dt = FT(30)
n = 50

# Specify the domain boundaries
zmax = FT(0)
zmin = FT(-1)
domain = Domains.IntervalDomain(zmin, zmax, x3boundary = (:bottom, :top))
mesh = Meshes.IntervalMesh(domain, nelems = n)

cs = Spaces.CenterFiniteDifferenceSpace(mesh)
fs = Spaces.FaceFiniteDifferenceSpace(cs)
zc = Fields.coordinate_field(cs)


#Parameter structure
p = [msp, param_set, zc]

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
Y = ArrayPartition(θ_l, θ_i, ρe_int)
function ∑tendencies!(dY, Y, p, t)
    compute_integrated_rhs!(dY, Y,t, p)
end
prob = ODEProblem(∑tendencies!, Y, (t0, tf),p)
sol = solve(
    prob,
    TsitPap8(),
    dt = dt,
    saveat = 60 * dt,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
);


#Plotting
plot1 = plot(
    parent(sol.u[1].x[1]),parent(zc),
    xlim = (0, ν),
    ylim = (-1, 0),
    label = "IC",
    lc = :black,
    lw = 2,
    legend = :outerright,
    xlabel = "θ(z)",
    ylabel = "z",
)
plot!(
    parent(sol.u[36].x[1]),parent(zc),
    label = "18h",
    lc = :red,
    lw = 2,
 
)
plot!(
    parent(sol.u[72].x[1]),parent(zc),
    label = "36h",
    lc = :blue,
    lw = 2,
 
)
plot!(
    parent(sol.u[108].x[1]),parent(zc),
    label = "54h",
    lc = :orange,
    lw = 2,
 
)

plot!(
    parent(sol.u[end].x[1]),parent(zc),
    label = "72h",
    lc = :green,
    lw = 2,
 
)


ρc_s = volumetric_heat_capacity.(sol.u[1].x[1], sol.u[1].x[2], ρc_ds, Ref(param_set))
temp = temperature_from_ρe_int.(sol.u[1].x[3], sol.u[1].x[2], ρc_s, Ref(param_set))
plot2 = plot(
    parent(temp),parent(zc),
   # xlim = (0, ν),
    ylim = (-1, 0),
    label = "IC",
    lc = :black,
    lw = 2,
    legend = :outerright,
    xlabel = "T(K)",
    ylabel = "z",
)
ρc_s = volumetric_heat_capacity.(sol.u[36].x[1], sol.u[36].x[2], ρc_ds, Ref(param_set))
temp = temperature_from_ρe_int.(sol.u[36].x[3], sol.u[36].x[2], ρc_s, Ref(param_set))
plot!(
    parent(temp),parent(zc),
    label = "18h",
    lc = :red,
    lw = 2,
 
)
ρc_s = volumetric_heat_capacity.(sol.u[72].x[1], sol.u[72].x[2], ρc_ds, Ref(param_set))
temp = temperature_from_ρe_int.(sol.u[72].x[3], sol.u[72].x[2], ρc_s, Ref(param_set))
plot!(
    parent(temp),parent(zc),
    label = "36h",
    lc = :blue,
    lw = 2,
 
)
ρc_s = volumetric_heat_capacity.(sol.u[108].x[1], sol.u[108].x[2], ρc_ds, Ref(param_set))
temp = temperature_from_ρe_int.(sol.u[108].x[3], sol.u[108].x[2], ρc_s, Ref(param_set))
plot!(
    parent(temp),parent(zc),
    label = "54h",
    lc = :orange,
    lw = 2,
 
)
ρc_s = volumetric_heat_capacity.(sol.u[end].x[1], sol.u[end].x[2], ρc_ds, Ref(param_set))
temp = temperature_from_ρe_int.(sol.u[end].x[3], sol.u[end].x[2], ρc_s, Ref(param_set))
plot!(
    parent(temp),parent(zc),
    label = "72h",
    lc = :green,
    lw = 2,
 
)

plot(plot1,plot2)
