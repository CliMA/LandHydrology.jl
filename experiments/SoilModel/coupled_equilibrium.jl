import ClimaCore: Fields


using CLIMAParameters
using CLIMAParameters.Planet: ρ_cloud_liq, ρ_cloud_ice, cp_l, cp_i, T_0, LH_f0
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using OrdinaryDiffEq: ODEProblem, solve, SSPRK33
using Plots
using UnPack
using LandHydrology
using LandHydrology.Domains
using LandHydrology.SoilInterface
using LandHydrology.SoilInterface.SoilHeatParameterizations
using LandHydrology.SoilInterface.SoilWaterParameterizations

const FT = Float64

# General soil composition
ν = FT(0.395);
ν_ss_quartz = FT(0.92)
ν_ss_minerals = FT(0.08)
ν_ss_om = FT(0.0)
ν_ss_gravel = FT(0.0);

#Water specific
Ksat = FT(4.42 / 3600 / 100) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(1.89)
vg_α = FT(7.5); # inverse meters
vg_m = FT(1) - FT(1) / vg_n
θ_r = FT(0)

#Heat specific
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

#collect all params
msp = SoilParams{FT}(
    ν,
    vg_n,
    vg_α,
    vg_m,
    Ksat,
    θ_r,
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
tf = FT(60 * 60 * 72)
dt = FT(30)
n = 50

zmax = FT(0)
zmin = FT(-1)
domain = Column(FT, zlim = (zmin, zmax), nelements = n)


#Boundary conditions
top_water_flux = FT(0)
top_heat_flux = FT(0)
bottom_water_flux = FT(0)
bottom_heat_flux = FT(0)
bc = SoilDomainBC(
    top = SoilComponentBC(
        hydrology = VerticalFlux(top_water_flux),
        energy = VerticalFlux(top_heat_flux),
    ),
    bottom = SoilComponentBC(
        hydrology = VerticalFlux(bottom_water_flux),
        energy = VerticalFlux(bottom_heat_flux),
    ),
)


# create model
soil_model = SoilModel(
    domain = domain,
    energy_model = SoilEnergyModel(),
    hydrology_model = SoilHydrologyModel(),
    boundary_conditions = bc,
    soil_param_set = msp,
    earth_param_set = param_set,
)

# initial conditions

function ic(z, t0, model)
    soil_param_set = model.soil_param_set
    earth_param_set = model.earth_param_set
    c = 20.0
    zmin = -1.0
    zmax = 0.0
    θ_max = 0.1975
    θ_min = 0.158
    θ_i = 0.0
    θ_l = θ_min + (θ_max - θ_min) * exp(-(z - zmax) / (zmin - zmax) * c)
    T = 288.0 + (289.0 - 288.0) * exp(-(z - zmax) / (zmin - zmax) * c)
    ρc_ds = soil_param_set.ρc_ds
    ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρc_ds, earth_param_set)
    ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, earth_param_set)
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
);


#Plotting
space_c, _ = make_function_space(domain)
zc = Fields.coordinate_field(space_c)
z = parent(zc)
ϑ_l = [parent(sol.u[k].ϑ_l) for k in 1:length(sol.u)]
θ_i = [parent(sol.u[k].θ_i) for k in 1:length(sol.u)]
ρe_int = [parent(sol.u[k].ρe_int) for k in 1:length(sol.u)]
#function expected(z, model)
#    @unpack ν, vgα, vgm, vgn  =   model.soil_param_set
#    C = -0.985118
#    return ν*(1.0+(vgα*abs(-z +C))^vgn)^(-vgm)
#end




indices = [1, 36, 72, 108, length(sol.t)]
labels = ["IC", "18h", "36h", "54h", "72h"]
plot1 = plot(
    xlim = (0.0, 0.4),
    ylim = (-1, 0),
    legend = :outerright,
    xlabel = "θ(z)",
    ylabel = "z",
)
for i in 1:1:length(indices)
    plot!(ϑ_l[indices[i]], z, label = labels[i], lw = 2)
end
#plot!(expected.(z,Ref(soil_model)),z, label= "Expected")
plot2 =
    plot(ylim = (-1, 0), legend = :outerright, xlabel = "T(K)", ylabel = "z")
for i in 1:1:length(indices)
    k = indices[i]
    ρc_s = volumetric_heat_capacity.(ϑ_l[k], θ_i[k], ρc_ds, Ref(param_set))
    temp = temperature_from_ρe_int.(ρe_int[k], θ_i[k], ρc_s, Ref(param_set))
    plot!(temp, z, label = labels[i], lw = 2)
end

plot(plot1, plot2)
