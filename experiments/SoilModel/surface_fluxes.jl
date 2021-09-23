import ClimaCore: Fields
using CLIMAParameters

struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
using UnPack
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33

using LandHydrology

using LandHydrology.Domains: Column, make_function_space
using LandHydrology.SoilInterface
using LandHydrology.SoilInterface.SoilWaterParameterizations
using LandHydrology.SoilInterface.SoilHeatParameterizations
using Plots
using SurfaceFluxes: DGScheme, surface_conditions
using CLIMAParameters.Planet:
    R_v, grav, ρ_cloud_liq, R_d, cp_v, LH_v0, T_0, cp_d
using Thermodynamics:
    Liquid,
    q_vap_saturation_generic,
    saturation_vapor_pressure,
    PhasePartition,
    cp_m

const FT = Float64
ν = FT(0.55)
Ksat = FT(1.31 / 100 / 3600 / 1000) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(1.68)
vg_α = FT(5.0) # inverse meters
θ_r = FT(0.084)
ν_ss_quartz = FT(0.4)
ν_ss_om = FT(0.0)
ν_ss_gravel = FT(0.0)
ρp = FT(1770 / (1.0 - ν))
κ_quartz = FT(7.7) # W/m/K
κ_minerals = FT(2.5) # W/m/K
κ_om = FT(0.25) # W/m/K
κ_liq = FT(0.57) # W/m/K
κ_ice = FT(2.29) # W/m/K
κ_solid = k_solid(ν_ss_om, ν_ss_quartz, κ_quartz, κ_minerals, κ_om)
κ_sat_frozen = ksat_frozen(κ_solid, ν, κ_ice)
κ_sat_unfrozen = ksat_unfrozen(κ_solid, ν, κ_liq)
ρc_ds = FT((1 - ν) * 1.926e06)

msp = SoilParams{FT}(
    ν = ν,
    S_s = S_s,
    ν_ss_gravel = ν_ss_gravel,
    ν_ss_om = ν_ss_om,
    ν_ss_quartz = ν_ss_quartz,
    ρc_ds = ρc_ds,
    κ_solid = κ_solid,
    κ_sat_unfrozen = κ_sat_unfrozen,
    κ_sat_frozen = κ_sat_frozen,
)
hm = vanGenuchten{FT}(; n = vg_n, α = vg_α, Ksat = Ksat, θr = θ_r)

# Domain
n = 10

zmax = FT(0)
zmin = FT(-0.55)
domain = Column(FT, zlim = (zmin, zmax), nelements = n)

# Boundary conditions
T_surf = 299.0
ρ_a_sfc = 1.17
q_atm = 0.015
surface_bc = PrescribedAtmosForcing{FT}(
    u_atm = 0.34,
    θ_atm = T_surf,
    z_atm = 0.05,
    θ_scale = T_surf,
    ρ_a_sfc = ρ_a_sfc,
    q_atm = q_atm,
)

bottom_flux = FT(0)
bc = SoilColumnBC(;
    top = surface_bc,
    bottom = SoilComponentBC(
        energy = VerticalFlux(bottom_flux),
        hydrology = VerticalFlux(bottom_flux),
    ),
)

# The model
soil_model = SoilModel(
    FT;
    domain = domain,
    energy_model = SoilEnergyModel(),
    hydrology_model = SoilHydrologyModel{FT}(hydraulic_model = hm),
    boundary_conditions = bc,
    soil_param_set = msp,
    earth_param_set = param_set,
)

# initial conditions
function initial_conditions(z::ft, model::SoilModel) where {ft <: AbstractFloat}
    param_set = model.earth_param_set
    hm = model.hydrology_model.hydraulic_model
    sp = model.soil_param_set
    T = 298.5
    @unpack ν, S_s = sp
    θ_i = 0.0
    θ_l = hydrostatic_profile(hm, z, -0.55, ν, S_s)
    ρcds = sp.ρc_ds
    ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρcds, param_set)
    ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, param_set)
    return (ϑ_l = θ_l, θ_i = θ_i, ρe_int = ρe_int)
end

t0 = FT(0)
tf = FT(3600 * 24 * 480)
dt = FT(160)
Y, Ya = initialize_states(soil_model, initial_conditions, t0)
soil_rhs! = make_rhs(soil_model)
prob = ODEProblem(soil_rhs!, Y, (t0, tf), Ya)



# solve simulation
sol = solve(prob, SSPRK33(), dt = dt, saveat = 3600 * 4)
space_c, _ = make_function_space(domain)
zc = Fields.coordinate_field(space_c).z
N = length(sol.t)
tm = [parent(sol.u[k].soil.ϑ_l)[end] for k in 1:N]
_T_ref = FT(T_0(param_set))
function f(u_k)
    θ_l = u_k.soil.ϑ_l
    θ_i = u_k.soil.θ_i
    ρeint = u_k.soil.ρe_int
    ν_eff = ν .- θ_i
    vlf = volumetric_liquid_fraction.(θ_l, ν_eff)
    ρcds = soil_model.soil_param_set.ρc_ds
    ρc_s = volumetric_heat_capacity.(vlf, θ_i, ρcds, Ref(param_set))
    temp = temperature_from_ρe_int.(ρeint, θ_i, ρc_s, Ref(param_set))

    Y = Fields.FieldVector(ϑ_l = θ_l, θ_i = θ_i, T = temp)
    bf = boundary_fluxes(
        Y,
        soil_model.boundary_conditions.top,
        :top,
        soil_model,
        space_c,
        0.0,
    )
    E = bf.fϑ_l
    F = bf.fρe_int
    Tsfc = parent(temp)[end]
    h_d = cp_d(param_set) * (Tsfc - _T_ref) + R_d(param_set) * _T_ref
    lh = FT(cp_v(param_set) * (Tsfc - _T_ref) + LH_v0(param_set))
    shf_t = F - E * (lh - h_d) * ρ_cloud_liq(param_set)
    lhf_t = lh * E * ρ_cloud_liq(param_set)
    return (; F, E, Tsfc, temp, θ_l, ρeint, shf_t, lhf_t)
end
nt = map(u_k -> f(u_k), sol.u)
en = getproperty.(nt, :F)
hyd = getproperty.(nt, :E)
tt = getproperty.(nt, :Tsfc)
tprof = getproperty.(nt, :temp)
θprof = getproperty.(nt, :θ_l)
ρeprof = getproperty.(nt, :ρeint)
shf = getproperty.(nt, :shf_t)
lhf = getproperty.(nt, :lhf_t)

indices = Array(1:1:6) * Int(round(length(sol.t) / 6))
plot1a = plot(
    parent(tprof[1]),
    parent(zc),
    xlabel = "T(K)",
    ylabel = "Depth (m)",
    label = "t = 0",
)
plot!(parent(tprof[indices[1]]), parent(zc), label = "t=20")
plot!(parent(tprof[indices[2]]), parent(zc), label = "t=40")
plot!(parent(tprof[indices[3]]), parent(zc), label = "t=60")
plot!(parent(tprof[indices[4]]), parent(zc), label = "t=80")
plot!(parent(tprof[indices[5]]), parent(zc), label = "t=100")
plot!(parent(tprof[indices[6]]), parent(zc), label = "t=120")
plot!(legend = :bottomleft)

plot1b = plot(
    parent(θprof[1]),
    parent(zc),
    xlabel = "θ",
    ylabel = "Depth (m)",
    label = "t = 0",
)
plot!(parent(θprof[indices[1]]), parent(zc), label = "t=20")
plot!(parent(θprof[indices[2]]), parent(zc), label = "t=40")
plot!(parent(θprof[indices[3]]), parent(zc), label = "t=60")
plot!(parent(θprof[indices[4]]), parent(zc), label = "t=80")
plot!(parent(θprof[indices[5]]), parent(zc), label = "t=100")
plot!(parent(θprof[indices[6]]), parent(zc), label = "t=120")
plot!(legend = :bottomright)


plot2 = plot(
    parent(ρeprof[1]),
    parent(zc),
    xlabel = "ρe_int",
    ylabel = "Depth (m)",
    label = "t = 0",
)
plot!(parent(ρeprof[indices[1]]), parent(zc), label = "t=20")
plot!(parent(ρeprof[indices[2]]), parent(zc), label = "t=40")
plot!(parent(ρeprof[indices[3]]), parent(zc), label = "t=60")
plot!(parent(ρeprof[indices[4]]), parent(zc), label = "t=80")
plot!(parent(ρeprof[indices[5]]), parent(zc), label = "t=100")
plot!(parent(ρeprof[indices[6]]), parent(zc), label = "t=120")
plot!(legend = :bottomright)
plot(plot1a, plot1b, plot2)
savefig("./profiles.png")

plot4 = plot(
    sol.t / 3600 / 24,
    en,
    label = "total heat flux",
    ylabel = "(W/m^2)",
    xlabel = "days",
)
plot!(sol.t / 3600 / 24, shf, label = "shf")
plot!(sol.t / 3600 / 24, lhf, label = "lhf")
plot!(legend = :topright)

plot5 = plot(
    sol.t / 3600 / 24,
    tt,
    label = "T",
    ylabel = "T_surf (K)",
    xlabel = "days",
)
plot!([0, sol.t[end] / 3600 / 24], [299.0, 299.0], label = "T_atm")
plot!(legend = :bottomright)
plot(plot4, plot5)
savefig("./heat_fluxes.png")

q_sat = q_vap_saturation_generic.(Ref(param_set), tt, ρ_a_sfc, Ref(Liquid()))
S_l_eff = min.(effective_saturation.(ν, tm, θ_r), 1.0)
ψ = matric_potential.(Ref(hm), S_l_eff)
g = grav(param_set)
Rv = R_v(param_set)
correction = exp.(g .* ψ ./ Rv ./ tt)
q_soil = q_sat .* correction
# E ∝ (q - q atm) # no soil resistance yet
# This potential rate is an odd definition to me. It is the potential rate
# assuming the same drag coefficient. However, if q_surf changes
# the drag coefficient will change because L_mo will change.
# With current SurfaceFluxes.jl, though, where moisture does not affect L_mo,
# this should be correct. 
E0 = hyd ./ (q_soil .- q_atm) .* (q_sat .- q_atm)



plot3 = plot(sol.t ./ 3600 / 24, tm, label = "ϑ_l")
plot!(xlabel = "time (days)", ylabel = "vwc at top")
plot!([0, sol.t[end] / 3600 / 24], [θ_r, θ_r], ls = :dash, label = "θᵣ")
plot!(legend = :bottomleft)

plot6 = plot(
    sol.t / 3600 / 24,
    hyd * 1000 * 3600 * 24,
    label = "E (w/o resistance)",
    ylabel = "E (mm/day)",
    xlabel = "days",
)
plot!(sol.t / 3600 / 24, E0 * 1000 * 3600 * 24, label = "E(potential)")

plot7 = plot(sol.t / 3600 / 24, q_sat, label = "q_sat(T)", xlabel = "days")
plot!(sol.t / 3600 / 24, q_soil, label = "q_soil(T,ψ)")
plot!([0, sol.t[end] / 3600 / 24], [q_atm, q_atm], label = "q_atm")
plot!(legend = :bottomright)
plot(plot3, plot6, plot7)
savefig("./moisture_fluxes.png")
