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
using OrdinaryDiffEq:
    ODEProblem,
    solve,
    SSPRK33,
    Rosenbrock23,
    Tsit5,
    SSPRK432,
    Feagin14,
    TsitPap8,
    CarpenterKennedy2N54
using DifferentialEquations
using Plots
using DelimitedFiles
using UnPack
using LandHydrology
using LandHydrology.SoilInterface
using LandHydrology.SoilHeatParameterizations
using LandHydrology.SoilWaterParameterizations
using Test
using Statistics

const FT = Float64

abstract type BC{FT <: AbstractFloat} end

struct FluxBC{FT} <: BC{FT}
    top_heat_flux::FT
    top_water_flux::FT
    btm_heat_flux::FT
    btm_water_flux::FT
end


function compute_integrated_rhs!(dY, Y, t, p)

    sp = p[1]
    param_set = p[2]
    zc = p[3]
    @unpack top_heat_flux, btm_heat_flux, top_water_flux, btm_water_flux = p[4]
    @unpack ν, vgn, vgα, vgm, ksat, θr, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp

    (Y_hydro, Y_energy) = Y.x
    (dY_hydro, dY_energy) = dY.x
    @unpack ϑ_l, θ_i = Y_hydro
    @unpack ρe_int = Y_energy

    dϑ_l = dY_hydro.ϑ_l
    dθ_i = dY_hydro.θ_i
    dρe_int = dY_energy.ρe_int

    # Compute center values of everything
    θ_l = ϑ_l
    ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
    T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))
    κ_dry = k_dry(param_set, sp)
    S_r = relative_saturation.(θ_l, θ_i, ν)
    kersten = kersten_number.(θ_i, S_r, Ref(sp))
    κ_sat =
        saturated_thermal_conductivity.(θ_l, θ_i, κ_sat_unfrozen, κ_sat_frozen)
    κ = thermal_conductivity.(κ_dry, kersten, κ_sat)
    ρe_int_l = volumetric_internal_energy_liq.(T, Ref(param_set))

    cs = axes(θ_i)


    S = effective_saturation.(θ_l; ν = ν, θr = θr)
    K = hydraulic_conductivity.(S; vgm = vgm, ksat = ksat)
    ψ = matric_potential.(S; vgn = vgn, vgα = vgα, vgm = vgm)
    h = ψ .+ zc

    interpc2f = Operators.InterpolateC2F()
    gradc2f_heat = Operators.GradientC2F()
    gradf2c_heat = Operators.GradientF2C(
        top = Operators.SetValue(top_heat_flux),
        bottom = Operators.SetValue(btm_heat_flux),
    )

    gradc2f_water = Operators.GradientC2F()
    gradf2c_water = Operators.GradientF2C(
        top = Operators.SetValue(top_water_flux),
        bottom = Operators.SetValue(btm_water_flux),
    )

    @. dϑ_l = -gradf2c_water(-interpc2f(K) * gradc2f_water(h)) #Richards equation
    @. dρe_int =
        -gradf2c_heat(
            -interpc2f(κ) * gradc2f_heat(T) -
            interpc2f(ρe_int_l * K) * gradc2f_water(h),
        )
    dθ_i = Fields.zeros(eltype(θ_i), cs)

    return dY

end

# General composition
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
domain = Domains.IntervalDomain(zmin, zmax, x3boundary = (:bottom, :top))
mesh = Meshes.IntervalMesh(domain, nelems = n)

cs = Spaces.CenterFiniteDifferenceSpace(mesh)
fs = Spaces.FaceFiniteDifferenceSpace(cs)
zc = Fields.coordinate_field(cs)

#Boundary conditions
top_water_flux = FT(0)
top_heat_flux = FT(0)
bottom_water_flux = FT(0)
bottom_heat_flux = FT(0)

flux_bc =
    FluxBC(top_heat_flux, top_water_flux, bottom_heat_flux, bottom_water_flux)

#Parameter structure
p = [msp, param_set, zc, flux_bc]

#initial conditions

# clean up parameters
function init_hydrology(z)
    θ_i = 0.0
    c = (20.0)
    zmax = 0.0
    zmin = -1.0
    theta_max = 0.1975
    theta_min = 0.158
    ϑ_l =
        theta_min +
        (theta_max - theta_min) * exp(-(z - zmax) / (zmin - zmax) * c)
    return (ϑ_l = ϑ_l, θ_i = θ_i)
end


function init_energy(z, soil_params, global_params)
    T_max = 289.0
    T_min = 288.0
    c = 20.0
    zmax = 0.0
    zmin = -1.0
    T = T_min + (T_max - T_min) * exp(-(z - zmax) / (zmin - zmax) * c)
    θ_i = 0.0
    theta_max = 0.1975
    theta_min = 0.158
    θ_l =
        theta_min +
        (theta_max - theta_min) * exp(-(z - zmax) / (zmin - zmax) * c)

    ρc_ds = soil_params.ρc_ds
    ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρc_ds, global_params)
    ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, global_params)
    return (ρe_int = ρe_int,)
end


hydrology_model = SoilHydrologyModel(init_hydrology, nothing)
energy_model = SoilEnergyModel(init_energy, nothing)
soil_model = SoilModel(energy_model, hydrology_model, msp, param_set)
Y = init_prognostic_vars(soil_model, cs)

function ∑tendencies!(dY, Y, p, t)
    #intermediate step to be added if needed
    compute_integrated_rhs!(dY, Y, t, p)
end

SubModel1 = ... 
#SoilModel = NoSoil(T_prescribed = T(x,y,z,t), \theta_i = f)])
SoilModel = FancySoil(args...; kwargs...)
CoupledModel = ...
main_model = MainModel(SubModel1, SoilModel, CoupledModel, ..., args...; kwargs...)

SubIC1 = sub_ic1 # functions f(x, y, z, params; kwargs...) => f(x, y, z, params; T_sfc = FieldReader)
SoilIC = soil_ic
CoupledIC = coupled_ic
ic = (SubIC1, SoilIC, ..., args...; kwargs...) # creates prognostic and aux components: Y = (Y_prog = (sub1 = Y1, soil = Y_soil, coupled = Y_coupled)

function make_rhs(main_model)
    rhs_sub1! = make_rhs(main_model.sub1) # rhs_sub1! updates dY.sub1 in place
    rhs_soil! = make_rhs(main_model.soil) # ""
    rhs_coupled! = make_rhs(main_model.coupled) # ""
    function rhs!(dY, Y, p, t)
        rhs_sub1!(dY.soil, Y, p, t)
        rhs_soil!(dY.soil, Y, p, t)
        rhs_coupled!(dY.soil, Y, p, t)
    end

    return rhs!
end

function make_rhs!(main_model)
    rhs_list = []
    for submodel in submodels:
        push!(rhs_list, make_rhs(submodel))
    end

    function rhs!(dY, Y, p, t)
        @unroll for (submodel, rhs!) in zip(submodels, rhs_list)
            rhs!(
                getproperty(dY, submodel),
                getproperty(Y, submodel),
                p,
                t
            )
        end
    end

    return rhs!
end

# Y = (Y1, SoilY, ...)
state = make_state(ic) # assembles a state, populates stuff needs to know that for SoilIC - dont add anything

ic = ModelIC(canopy = canopy_ic, soil = soil_ic, river = river_ic,...)
sim = Simulation(model, ic, stepper, dt, ...)

prob = ODEProblem(∑tendencies!, Y, (t0, tf), p)

# solve simulation
sol = solve(
    prob,
    CarpenterKennedy2N54(),
    dt = dt,
    saveat = 60 * dt,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
);


#Plotting
z = parent(zc)
ϑ_l = [parent(sol.u[k].x[1].ϑ_l) for k in 1:length(sol.u)]
θ_i = [parent(sol.u[k].x[1].θ_i) for k in 1:length(sol.u)]
ρe_int = [parent(sol.u[k].x[2].ρe_int) for k in 1:length(sol.u)]
indices = [1, 36, 72, 108, length(sol.t)]
labels = ["IC", "18h", "36h", "54h", "72h"]
plot1 = plot(
    xlim = (0.1, 0.25),
    ylim = (-1, 0),
    legend = :outerright,
    xlabel = "θ(z)",
    ylabel = "z",
)
for i in 1:1:length(indices)
    plot!(ϑ_l[indices[i]], z, label = labels[i], lw = 2)
end
plot2 =
    plot(ylim = (-1, 0), legend = :outerright, xlabel = "T(K)", ylabel = "z")
for i in 1:1:length(indices)
    k = indices[i]
    ρc_s = volumetric_heat_capacity.(ϑ_l[k], θ_i[k], ρc_ds, Ref(param_set))
    temp = temperature_from_ρe_int.(ρe_int[k], θ_i[k], ρc_s, Ref(param_set))
    plot!(temp, z, label = labels[i], lw = 2)
end

plot(plot1, plot2)
