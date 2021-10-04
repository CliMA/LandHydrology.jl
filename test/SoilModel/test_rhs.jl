#  K = 0
# T < T_freeze in part of domain
# test that rhs = source
using ClimaCore: Fields
using CLIMAParameters
using CLIMAParameters.Planet: ρ_cloud_liq, T_freeze, grav, ρ_cloud_ice, LH_f0
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
t0 = FT(0)
tf = FT(60 * 60 * 24 * 32)
dt = FT(20)
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
Y = set_initial_state(soil_model, initial_conditions, 0.0)
soil_rhs! = make_rhs(soil_model)
dY = similar(Y) .+ 0.0
soil_rhs!(dY,Y,nothing,0.0)
ρ_liq = ρ_cloud_liq(param_set)
ρ_ice = ρ_cloud_ice(param_set)
@test sum(parent(dY_updated.ϑ_l .*ρ_liq .+ dY_updated.θ_i .*ρ_ice) .== 0.0) .== 20

computed_source! = add_source(PhaseChange{FT}(Δz = Δz), soil_model)
dY2 = similar(Y) .+ 0.0
computed_source!(dY2, Y, nothing, 0.0)
@test parent(dY2.ϑ_l) ≈ parent(dY.ϑ_l)
@test parent(dY2.θ_i) ≈ parent(dY.θ_i)
