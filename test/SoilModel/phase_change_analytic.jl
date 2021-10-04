using ClimaCore: Fields
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using SpecialFunctions
using UnPack
using SciMLBase: init, step!
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
using CLIMAParameters.Planet: ρ_cloud_liq, T_freeze, grav, ρ_cloud_ice, LH_f0
using Plots
using ArtifactWrappers
using DelimitedFiles
FT = Float64
ν = FT(0.535)
Ksat = FT(0.0)
S_s = FT(1e-3) #inverse meters
vg_n = FT(1.48)
vg_α = FT(1.11) # inverse meters
θ_r = FT(0.0)
ν_ss_quartz = FT(0.2)
ν_ss_minerals = FT(0.6)
ν_ss_om = FT(0.2)
ν_ss_gravel = FT(0.0)
κ_quartz = FT(7.7)
κ_minerals = FT(2.5)
κ_om = FT(0.25)
κ_liq = FT(0.57)
κ_ice = FT(2.29);
κ_solid = k_solid(ν_ss_om, ν_ss_quartz, κ_quartz, κ_minerals, κ_om)
κ_sat_frozen = ksat_frozen(κ_solid, ν, κ_ice)
κ_sat_unfrozen = ksat_unfrozen(κ_solid, ν, κ_liq);
ρp = FT(2700)
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
tf = FT(3600 * 20 *24)
dt = FT(50)
n = 40

zmax = FT(0)
zmin = FT(-3)
Δz = FT((zmax-zmin)/n)
domain = Column(FT, zlim = (zmin, zmax), nelements = n)

#Boundary conditions
top_flux = FT(0)
bottom_flux = FT(0)
bc = SoilDomainBC(
    domain;
    top = SoilComponentBC(
        hydrology = VerticalFlux(top_flux),
        energy = Dirichlet((t)-> eltype(t)(263.15))
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
    T = 275.15
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

#integrator = init(prob, SSPRK33(); dt = dt)
function get_T(vlf, vif, ρeint)
    ρcs = volumetric_heat_capacity.(vlf, vif, ρc_ds, Ref(param_set))
    
    temp = temperature_from_ρe_int.(ρeint, vif, ρcs, Ref(param_set)) .- 273.15
    return temp
end

function make_plot(integrator)
    iter = 1:1:length(integrator.sol.t)
    ρe = [parent(integrator.sol[k].ρe_int)[end] for k in iter]
    vlf = [parent(integrator.sol[k].ϑ_l)[end] for k in iter]
    vif = [parent(integrator.sol[k].θ_i)[end] for k in iter]
    T = get_T.(vlf, vif, ρe)
    plot1= plot(T)
    plot2 = plot(vif, label = "ice")
    plot!(vlf, label = "liq")
    plot(plot1, plot2)
    
end    
space_c, _ = make_function_space(domain)
zc = Fields.coordinate_field(space_c)
z = parent(zc)[:]

# solve simulation
sol = solve(
    prob,
    SSPRK33(),
    dt = dt,
    saveat = 3600,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
)


# Compute the thermal conductivity and heat capacity in the frozen region - subscript 1.
θ_l0 = FT(0.33)
kdry = k_dry(param_set, msp)
ksat = saturated_thermal_conductivity(
    FT(0.0),
    θ_l0 * ρ_cloud_liq(param_set) / ρ_cloud_ice(param_set),
    κ_sat_unfrozen,
    κ_sat_frozen,
)
kersten = kersten_number(
    θ_l0 * ρ_cloud_liq(param_set) / ρ_cloud_ice(param_set),
    θ_l0 * ρ_cloud_liq(param_set) / ρ_cloud_ice(param_set) / ν,
    msp,
)
λ1 = thermal_conductivity(kdry, kersten, ksat)
c1 = volumetric_heat_capacity(
    FT(0.0),
    θ_l0 * ρ_cloud_liq(param_set) / ρ_cloud_ice(param_set),
    ρc_ds,
    param_set,
)
d1 = λ1 / c1;

# Compute the thermal conductivity and heat capacity in the region
# with liquid water - subscript 2.
ksat =
    saturated_thermal_conductivity(θ_l0, FT(0.0), κ_sat_unfrozen, κ_sat_frozen)
kersten = kersten_number(FT(0.0), θ_l0 / ν, msp)
λ2 = thermal_conductivity(kdry, kersten, ksat)
c2 = volumetric_heat_capacity(θ_l0, FT(0.0), ρc_ds, param_set)
d2 = λ2 / c2;

# Initial T and surface T, in Celsius
Ti = FT(2)
Ts = FT(-10.0);

# The solution requires the root of the implicit equation below
# (you'll need to install `Roots` via the package manager).
# ``` julia
# using Roots
# function implicit(ζ)
#    term1 = exp(-ζ^2) / ζ / erf(ζ)
#    term2 =
#        -λ2 * sqrt(d1) * (Ti - 0) /
#        (λ1 * sqrt(d2) * (0 - Ts) * ζ * erfc(ζ * sqrt(d1 / d2))) *
#        exp(-d1 / d2 * ζ^2)
#    term3 =
#        -LH_f0(param_set) * ρ_cloud_liq(param_set) * θ_l0 * sqrt(π) / c1 /
#        (0 - Ts)
#    return (term1 + term2 + term3)
# end
# find_zero(implicit, (0.25,0.27), Bisection())
# ```

# The root is
ζ = 0.26447353269809687;



# This function plots the analytic solution
# and the simulated result, at the output time index
# `k`.
function f(k; ζ = ζ)
    vlf = parent(sol.u[k].ϑ_l)[:]
    ρeint = parent(sol.u[k].ρe_int)[:]
    vif = parent(sol.u[k].θ_i)[:]
    ρcs = volumetric_heat_capacity.(vlf, vif, ρc_ds, Ref(param_set))
    
    temp = temperature_from_ρe_int.(ρeint, vif, ρcs, Ref(param_set)) .- 273.15
    plot(
        temp,
        z,
        xlim = [-10.5, 3],
        ylim = [zmin, zmax],
        xlabel = "T",
        label = "simulation",
    )
    t = sol.t[k]
    zf = 2.0 * ζ * sqrt(d1 * t)
    myz = zmin:0.001:0
    spatially_varying =
        (erfc.(abs.(myz) ./ (zf / ζ / (d1 / d2)^0.5))) ./
        erfc(ζ * (d1 / d2)^0.5)
    mask = abs.(myz) .>= zf
   
    plot!(
        Ti .- (Ti - 0.0) .* spatially_varying[mask],
        myz[mask],
        label = "analytic",
        color = "green",
    )
    spatially_varying = ((erf.(abs.(myz) ./ (zf / ζ)))) ./ erf(ζ)
    mask = abs.(myz) .< zf
    plot!(
        Ts .+ (0.0 - Ts) .* spatially_varying[mask],
        myz[mask],
        label = "",
        color = "green",
    )
end

# Now we will plot this and compare to other methods for modeling phase change without water movement.
# These solutions were produced by modifying the Supplemental Program 5:3 from [Bonan19a](@cite).

# Excess heat solution:
eh_dataset = ArtifactWrapper(
    @__DIR__,
    isempty(get(ENV, "CI", "")),
    "eh",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/6xs1r98wk7u1b0xjhpkdvt80q4sh16un.csv",
        filename = "bonan_data.csv",
    ),],
);
eh_dataset_path = get_data_folder(eh_dataset);
eh_data = joinpath(eh_dataset_path, "bonan_data.csv")
ds = readdlm(eh_data, ',');
# Apparent heat capacity solution:
ahc_dataset = ArtifactWrapper(
    @__DIR__,
    isempty(get(ENV, "CI", "")),
    "ahc",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/d6xciskzl2djwi03xxfo8wu4obrptjmy.csv",
        filename = "bonan_data_ahc.csv",
    ),],
);
ahc_dataset_path = get_data_folder(ahc_dataset);
ahc_data = joinpath(ahc_dataset_path, "bonan_data_ahc.csv")
ds2 = readdlm(ahc_data, ',');

k = length(sol.t)
f(k)
plot!(ds[:, 2], ds[:, 1], label = "Excess Heat")
plot!(ds2[:, 2], ds2[:, 1], label = "Apparent heat capacity")
plot!(legend = :bottomleft)
savefig("analytic_comparison.png")

