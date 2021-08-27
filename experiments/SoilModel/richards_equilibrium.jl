import ClimaCore:
    Fields,
    Domains,
    Meshes,
    Spaces

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

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
using LandHydrology.SoilWaterParameterizations
using Test
using Statistics

const FT = Float64

# General soil composition
ν = FT(0.495);

#Water specific
Ksat = FT(0.0443 / 3600 / 100) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(2)
vg_α = FT(2.6); # inverse meters
vg_m = FT(1) - FT(1) / vg_n
θ_r = FT(0)

#collect all params
msp = SoilWaterParams{FT}(
    ν,
    vg_n,
    vg_α,
    vg_m,
    Ksat,
    θ_r,
    S_s,
)


#Simulation and domain info
t0 = FT(0)
tf = FT(60 * 60 * 24 *36)
dt = FT(100)
n = 50

zmax = FT(0)
zmin = FT(-10)
domain = Domains.IntervalDomain(zmin, zmax, x3boundary = (:bottom, :top))
mesh = Meshes.IntervalMesh(domain, nelems = n)

cs = Spaces.CenterFiniteDifferenceSpace(mesh)
fs = Spaces.FaceFiniteDifferenceSpace(cs)
zc = Fields.coordinate_field(cs)

#Boundary conditions
top_water_flux = FT(0)
bottom_water_flux = FT(0)
water_bc = FluxBC(top_water_flux,bottom_water_flux)
bc = SoilBC(hydrology_bc = water_bc)

# create model
soil_model = RichardsEquation(bc, msp, param_set)

# initial conditions
function hydrology_ic(z, model)
    θ_i = 0.0
    θ_l  = 0.494
    return (ϑ_l = θ_l, θ_i = θ_i)
end

ic = SoilIC(hydrology = IC(hydrology_ic))

Y = create_initial_state(soil_model, ic, zc)
soil_rhs! = make_rhs(soil_model)
dY =  similar(Y)
soil_rhs!(dY,Y,[],0.0)
prob = ODEProblem(soil_rhs!, Y, (t0, tf), [])

# solve simulation
sol = solve(
    prob,
    #TsitPap8(), doesnt work, ArgumentError: broadcasting over dictionaries and `NamedTuple`s is reserved
    CarpenterKennedy2N54(), # this does work
    dt = dt,
    saveat = 60 * dt,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
);


#Plotting
z = parent(zc)
ϑ_l = [parent(sol.u[k].x[1].ϑ_l) for k in 1:length(sol.u)]

indices = [1, 88, length(sol.t)]
labels = ["IC", "6d", "36d"]
plot1 = plot(
    xlim = (0.4, 0.525),
    ylim = (-10, 0),
    legend = :outerright,
    xlabel = "θ(z)",
    ylabel = "z",
)
for i in 1:1:length(indices)
    plot!(ϑ_l[indices[i]], z, label = labels[i], lw = 2)
end



function expected(z, z_interface)
    ν = 0.495
    S_s = 1e-3
    α = 2.6
    n = 2.0
    m = 0.5
    if z < z_interface
        return -S_s * (z - z_interface) + ν
    else
        return ν * (1 + (α * (z - z_interface))^n)^(-m)
    end
end
plot!(expected.(z, -0.56), z, lw =1, label = "expected");
plot(plot1)
#put expected on plot
