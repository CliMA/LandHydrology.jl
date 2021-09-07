using ClimaCore: Fields
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using OrdinaryDiffEq:
    ODEProblem,
    solve,
    CarpenterKennedy2N54,# does not work
    SSPRK33
using LandHydrology
using LandHydrology.Domains
using LandHydrology.SoilInterface
using LandHydrology.SoilInterface.SoilWaterParameterizations
using LandHydrology.SoilInterface.SoilHeatParameterizations
using Statistics
using Test
@testset "Variably saturated equilibrium" begin
    FT = Float64

    # General soil composition
    ν = FT(0.495);
    #Water specific
    Ksat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6); # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0)
    
    ρc_ds = FT(2e6)
    #collect all params
    msp = SoilParams{FT}(
        ν,
        vg_n,
        vg_α,
        vg_m,
        Ksat,
        θ_r,
        S_s,
        0.0,
        0.0,
        0.0,
        ρc_ds,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0
    )
    
    
    #Simulation and domain info
    t0 = FT(0)
    tf = FT(60 * 60 * 24 *36)
    dt = FT(100)
    n = 50
    
    zmax = FT(0)
    zmin = FT(-10)
    domain = Column(FT, zlim = (zmin, zmax), nelements = n)
    
    #Boundary conditions
    top_water_flux = FT(0)
    bottom_water_flux = FT(0)
    bc = SoilDomainBC(top = SoilComponentBC(hydrology = VerticalFlux(top_water_flux)),
                      bottom = SoilComponentBC(hydrology = VerticalFlux(bottom_water_flux)))
    
    # create model
    soil_model = SoilModel(domain = domain,energy_model = PrescribedTemperatureModel(), hydrology_model = SoilHydrologyModel(),
                       boundary_conditions = bc, soil_param_set = msp, earth_param_set = param_set)
    
    # initial conditions
    function initial_conditions(z, t0, model)
        T = model.energy_model.T_profile(z,t0) # to be consistent with PrescribedT Default. 
        θ_i = 0.0
        θ_l = 0.494
        ρc_ds =model.soil_param_set.ρc_ds
        ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρc_ds, model.earth_param_set)
        ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, model.earth_param_set)
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
        saveat = 60 * dt,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
);
    
    space_c, _ = make_function_space(domain)
    zc = Fields.coordinate_field(space_c)
    z = parent(zc)
    ϑ_l = [parent(sol.u[k].ϑ_l) for k in 1:length(sol.u)]
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
    
    
    
    @test sqrt(mean(ϑ_l[end] .- expected.(z, -0.56)).^2.0) < 1e-4
end
#= 
#Plotting
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

plot!(expected.(z, -0.56), z, lw =1, label = "expected");
plot(plot1)
=#
