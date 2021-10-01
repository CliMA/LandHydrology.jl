using ClimaCore: Fields
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using UnPack
using OrdinaryDiffEq:
    ODEProblem,
    solve,
    CarpenterKennedy2N54,# does not work
    SSPRK33,
    SSPRK73
using LandHydrology
using LandHydrology.Domains: Column, make_function_space
using LandHydrology.SoilInterface
using LandHydrology.SoilInterface.SoilWaterParameterizations
using LandHydrology.SoilInterface.SoilHeatParameterizations

using ArtifactWrappers
using DelimitedFiles
FT = Float64
ν = FT(0.535)
Ksat = FT(3.2e-6)
S_s = FT(1e-3) #inverse meters
vg_n = FT(1.48)
vg_α = FT(1.11) # inverse meters
θ_r = FT(0.05)
ν_ss_quartz = FT(0.7)
ν_ss_om = FT(0.3)
ν_ss_gravel = FT(0.0)

κ_quartz = FT(7.7)
κ_minerals = FT(2.4)
κ_om = FT(0.25)
κ_liq = FT(0.57)
κ_ice = FT(2.29);
κ_solid = k_solid(ν_ss_om, ν_ss_quartz, κ_quartz, κ_minerals, κ_om)
κ_sat_frozen = ksat_frozen(κ_solid, ν, κ_ice)
κ_sat_unfrozen = ksat_unfrozen(κ_solid, ν, κ_liq);
ρp = FT(3200)
ρc_ds = FT((1 - ν) * 2.3e6);

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
tf = FT(3600)#*50)
dt = FT(1)
n = 10

zmax = FT(0)
zmin = FT(-0.2)
Δz = FT((zmax-zmin)/n)
domain = Column(FT, zlim = (zmin, zmax), nelements = n)

#Boundary conditions
top_flux = FT(0)
bottom_flux = FT(0)
bc = SoilDomainBC(
    domain;
    top = SoilComponentBC(
        hydrology = VerticalFlux(top_flux),
        energy = ConductiveSurface{FT}(28.0, 267.15),
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
    sources = PhaseChange{FT}(Δz = Δz, γ_LTE = 10.0),
    soil_param_set = msp,
    earth_param_set = param_set,
)

# initial conditions
function initial_conditions(z::FT, t0::FT, model::SoilModel)
    param_set = model.earth_param_set
    T = 279.85
    θ_i = 0.0
    θ_l = 0.33
    ρcds = model.soil_param_set.ρc_ds
    ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρcds, param_set)
    ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, param_set)
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
    saveat = 3600,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
)

space_c, _ = make_function_space(domain)
zc = Fields.coordinate_field(space_c)
z = parent(zc)
vlf = parent(sol.u[end].ϑ_l)
ρeint = parent(sol.u[end].ρe_int)
ρc_s = volumetric_heat_capacity.(vlf, 0.0, ρc_ds, Ref(param_set))
#=
temp = temperature_from_ρe_int.(ρeint, 0.0, ρc_s, Ref(param_set))
dataset = ArtifactWrapper(
    @__DIR__,
    isempty(get(ENV, "CI", "")),
    "mizoguchi",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/3xbo4rlam8u390vmucc498cao6wmqlnd.csv",
        filename = "mizoguchi_all_data.csv",
    ),],
);
dataset_path = get_data_folder(dataset);
data = joinpath(dataset_path, "mizoguchi_all_data.csv")
ds = readdlm(data, ',')
hours = ds[:, 1][2:end]
vwc = ds[:, 2][2:end] ./ 100.0
depth = ds[:, 3][2:end]
mask_12h = hours .== 12
mask_24h = hours .== 24
mask_50h = hours .== 50;

plot_12h =
    scatter(vwc[mask_12h], -depth[mask_12h], label = "", color = "purple")
plot!(
    dons_arr[13]["soil.water.θ_i"] .+ dons_arr[13]["soil.water.ϑ_l"],
    z,
    label = "",
    lc = :green,
    lw = 2,
)
plot!(
    if_dons_arr[13]["soil.water.θ_i"] .+ if_dons_arr[13]["soil.water.ϑ_l"],
    z,
    label = "",
    lc = :orange,
    lw = 2,
)

plot!(title = "12h")
plot!(xlim = [0.2, 0.55])
plot!(xticks = [0.2, 0.3, 0.4, 0.5])
plot!(ylabel = "Depth (m)");

plot_24h =
    scatter(vwc[mask_24h], -depth[mask_24h], label = "Data", color = "purple")
plot!(
    dons_arr[25]["soil.water.θ_i"] .+ dons_arr[25]["soil.water.ϑ_l"],
    z,
    label = "Sim, Ω = 7",
    lc = :green,
    lw = 2,
)

plot!(
    if_dons_arr[25]["soil.water.θ_i"] .+ if_dons_arr[25]["soil.water.ϑ_l"],
    z,
    label = "Sim, Ω = 0",
    lc = :orange,
    lw = 2,
)

plot!(title = "24h")
plot!(legend = :bottomright)
plot!(xlim = [0.2, 0.55])
plot!(xticks = [0.2, 0.3, 0.4, 0.5]);

plot_50h =
    scatter(vwc[mask_50h], -depth[mask_50h], label = "", color = "purple")
plot!(
    dons_arr[51]["soil.water.θ_i"] .+ dons_arr[51]["soil.water.ϑ_l"],
    z,
    label = "",
    lc = :green,
    lw = 2,
)

plot!(
    if_dons_arr[51]["soil.water.θ_i"] .+ if_dons_arr[51]["soil.water.ϑ_l"],
    z,
    label = "",
    lc = :orange,
    lw = 2,
)

plot!(title = "50h")
plot!(xlim = [0.2, 0.55])
plot!(xticks = [0.2, 0.3, 0.4, 0.5]);

plot(plot_12h, plot_24h, plot_50h, layout = (1, 3))
plot!(xlabel = "θ_l+θ_i")
=#
