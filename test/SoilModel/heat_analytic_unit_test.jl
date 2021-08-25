# Test heat equation agrees with analytic solution to problem 55
# on page 28 in
# https://ocw.mit.edu/courses/mathematics/18-303-linear-partial-differential-equations-fall-2006/lecture-notes/heateqni.pdf.
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
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using OrdinaryDiffEq:
    ODEProblem,
    solve,
    CarpenterKennedy2N54
using DelimitedFiles
using UnPack
using LandHydrology
using LandHydrology.SoilInterface
using LandHydrology.SoilHeatParameterizations
using Test
using Statistics

const FT = Float64

abstract type bc end

struct TDirichlet{FiT} <: bc
    T::FiT
end

function compute_soil_heat_rhs!(
    dY,
    Y,
    t,
    p,
    top::TDirichlet,
    bottom::TDirichlet,
)
    sp = p[1]
    param_set = p[2]
    @unpack ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp

    (Y_energy,) = Y.x
    ρe_int = Y_energy.ρe_int
    (dY_energy,) = dY.x
    dρe_int = dY_energy.ρe_int
    
    θ_l, θ_i = 0.0, 0.0
    
    ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρc_ds, param_set)
    Ttop = top.T(t)
    Tbot = bottom.T(t)
    T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))

    κ_dry = k_dry(param_set, sp)
    S_r = relative_saturation(θ_l, θ_i, ν)
    kersten = kersten_number(θ_i, S_r, sp)
    κ_sat =
        saturated_thermal_conductivity(θ_l, θ_i, κ_sat_unfrozen, κ_sat_frozen)
    κ = thermal_conductivity(κ_dry, kersten, κ_sat)

    gradc2f = Operators.GradientC2F(
        top = Operators.SetValue(Ttop),
        bottom = Operators.SetValue(Tbot),
    )
    gradf2c = Operators.GradientF2C()

    @. dρe_int = gradf2c(κ * gradc2f(T)) # κ is just a constant here, eventually could be a field
    return dY

end
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
msp = SoilHeatParams(
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
    κ_dry_parameter,
)


t0 = FT(0)
tf = FT(2)
dt = FT(1e-4)
n = 60

# Specify the domain boundaries
zmax = FT(1)
zmin = FT(0)
domain = Domains.IntervalDomain(zmin, zmax, x3boundary = (:bottom, :top))
mesh = Meshes.IntervalMesh(domain, nelems = n)

cs = Spaces.CenterFiniteDifferenceSpace(mesh)
fs = Spaces.FaceFiniteDifferenceSpace(cs)
zc = Fields.coordinate_field(cs)
function init_energy(z, soil_params, global_params)
    T0 = 0.0
    θ_l, θ_i = 0.0, 0.0
    ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρc_ds, param_set)
    ρe_int = volumetric_internal_energy.(θ_i, ρc_s, T0, Ref(param_set))
    return (ρe_int = ρe_int,)
end

energy_model = SoilEnergyModel(init_energy, nothing)
hydrology_model = PrescribedHydrologyModel()
soil_model  = SoilModel(energy_model, hydrology_model, msp, param_set)
Y = init_prognostic_vars(soil_model, cs)
tau = FT(1) # period (sec)
A = FT(5) # amplitude (K)
ω = FT(2 * pi / tau)
topbc = TDirichlet(t -> eltype(t)(0.0))
bottombc = TDirichlet(t -> A * cos(ω * t))

p = [msp, param_set, topbc, bottombc]
function ∑tendencies!(dY, Y, p, t)
    top = p[3]
    bot = p[4]
    compute_soil_heat_rhs!(dY,Y, t, p, top, bot)
end

prob = ODEProblem(∑tendencies!, Y, (t0, tf), p)
sol = solve(
    prob,
    CarpenterKennedy2N54(),
    dt = dt,
)

t = sol.t
z = parent(zc)[:]
num =
    exp.(sqrt(ω / 2) * (1 + im) * (1 .- z)) .-
    exp.(-sqrt(ω / 2) * (1 + im) * (1 .- z))
denom = exp(sqrt(ω / 2) * (1 + im)) - exp.(-sqrt(ω / 2) * (1 + im))
analytic_soln = real(num .* A * exp(im * ω * tf) / denom)

ρe_int_final = parent(sol.u[end].x[1].ρe_int)
θ_l, θ_i = 0.0, 0.0
ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρc_ds, param_set)
Tfinal = temperature_from_ρe_int.(ρe_int_final, θ_i, ρc_s, Ref(param_set))
MSE = mean((analytic_soln .- Tfinal) .^ 2.0)
@test MSE < 1e-2
