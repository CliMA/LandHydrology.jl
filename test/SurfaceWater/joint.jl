using LandHydrology: PrescribedAtmosState
using LandHydrology.SurfaceWater: SurfaceWaterModel
FT = Float64

# General soil composition
ν = FT(0.495)
#Water specific
Ksat = FT(0.0443 / 3600 / 100) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(2.0)
vg_α = FT(2.6) # inverse meters
vg_m = FT(1) - FT(1) / vg_n
θ_r = FT(0)

#collect all params
msp = SoilParams{FT}(ν = ν, S_s = S_s)


#Simulation and domain info
t0 = FT(0)
tf = FT(60 * 60 * 24 * 36)
dt = FT(100)
n = 50

zmax = FT(0)
zmin = FT(-10)
domain = Column(FT, zlim = (zmin, zmax), nelements = n)

#Boundary conditions
top_water_flux = FT(0)
bottom_water_flux = FT(0)
bc = SoilColumnBC(;
                  top = SoilComponentBC(hydrology = VerticalFlux(top_water_flux)),
                  bottom = SoilComponentBC(hydrology = VerticalFlux(bottom_water_flux)),
                  )

# create model
hydraulics_model =
    vanGenuchten{FT}(n = vg_n, α = vg_α, Ksat = Ksat, θr = θ_r)

soil_model = SoilModel(
    FT;
    domain = domain,
    energy_model = PrescribedTemperatureModel(),
    hydrology_model = SoilHydrologyModel{FT}(
        hydraulic_model = hydraulics_model,
    ),
    boundary_conditions = bc,
    soil_param_set = msp,
    earth_param_set = param_set,
)
surface = SurfaceWaterModel{FT}()
function initial_conditions_sfc(z, model)
    return (;h = 0.0)
end

precip = FT(-1e-10)
land_model = LandHydrologyModel{FT}(soil_model,surface; atmos_state = PrescribedAtmosState{FT}(precip))
# initial conditions
@test_throws ErrorException default_initial_conditions(soil_model)
@test_throws ErrorException default_initial_conditions(land_model)
function initial_conditions(z, model)
    θ_i = 0.0
    θ_l = 0.494
    return (ϑ_l = θ_l, θ_i = θ_i)
end
Y, Ya = initialize_states(land_model, (;soil =initial_conditions,sfc_water = initial_conditions_sfc))
land_rhs! = make_rhs(land_model)
land_sim = Simulation(
    land_model,
    SSPRK33(),
    Y_init = Y,
    dt = dt,
    tspan = (t0, tf),
    Ya_init = Ya,
    saveat = 60 * dt,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
)

# solve simulation
@test step!(land_sim) isa Nothing # either error or integration runs
run!(land_sim)
sol = land_sim.integrator.sol

z = parent(Ya.soil.zc)

