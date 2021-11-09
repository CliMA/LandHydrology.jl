import XLSX
using OrdinaryDiffEq
using DifferentialEquations
using Plots
using SymPy
using LsqFit
using NLsolve
using Revise

xf = XLSX.readxlsx("/Users/annagl/Documents/LandHydrology.jl/src/parameters_roots.xlsx")
XLSX.sheetnames(xf)
sh = xf["Sheet1"] 

########## Parameters ###########
# need TOML pkg
#In parameter_file.toml
#[parameter_1]
#RunValue = 0.3
#Type = "real"

#[parameter_2]
#RunValue = 0.3
#Type = "real"
#In your julia file do_some_canopy.jl:
#using TOML
#filename = "parameter_file.toml"
#toml_dict = TOML.parsefile(filename)
#edited)
#5:28
#access via toml_dict["parameter_1"]["RunValue"] for example

MPa_to_Pa = sh[2,3]
rho_water = sh[3,3]
rhog_MPa = sh[5,3]
mass_mole_water = sh[6,3]
volume_mole_water = sh[7,3]
h_root = sh[8,3]
h_stem = sh[9,3]
LMA = sh[11,3]
WD = sh[12,3]
K_max_stem = sh[38,3]
K_max_root = sh[39,3]
K_max_stem_moles = sh[43,3]
K_max_root_moles = sh[44,3]
K_max_root_total_moles = sh[45,3]
porosity_leaf = sh[46,3]
porosity_stem = sh[47,3]
size_reservoir_stem = sh[50,3]
size_reservoir_leaf = sh[52,3]
size_reservoir_stem_moles = sh[54,3]
size_reservoir_leaf_moles = sh[55,3]
C_stem = sh[56,3]
C_leaf = sh[58,3]
S_stem = sh[59,3]
S_leaf = sh[60,3]
epsilon_leaf = sh[61,3] 
epsilon_stem = sh[62,3]
pi_0_leaf = sh[63,3]
pi_0_stem = sh[64,3]
pi_tlp_leaf = sh[65,3]
pi_tlp_stem = sh[66,3]
theta_r_leaf = sh[67,3]
theta_r_stem = sh[68,3]
theta_tlp_leaf = sh[69,3]
theta_tlp_stem = sh[70,3]
theta_ft_stem = sh[71,3]
theta_ft_leaf = sh[72,3]
f_cap = sh[73,3]
P_50_stem = sh[74,3] 
P_50_root = sh[76,3]
a_x_stem = sh[75,3]
a_x_root = sh[77,3]
a_root = sh[78,3]
a_stem = sh[79,3] 
b_root = sh[80,3]
b_stem = sh[81,3]

paramset = [a_root a_stem b_root b_stem K_max_root_total_moles K_max_stem_moles rho_water rhog_MPa h_root h_stem pi_0_leaf pi_0_stem epsilon_leaf epsilon_stem theta_r_leaf theta_r_stem WD theta_tlp_leaf theta_tlp_stem theta_ft_stem size_reservoir_stem_moles size_reservoir_leaf_moles theta_ft_leaf volume_mole_water S_stem S_leaf porosity_stem porosity_leaf];

# Functions
 
function vc_integral_approx(z_do, z_up, p_do, p_up, a, b, Kmax)
    num_do = log.(a .* exp.(b.*p_do) .+ 1);
    num_up = log.(a .* exp.(b.*p_up) .+ 1);
    flow_approx = -Kmax .* (a+1) / (a*b) .* (num_up .- num_do) .* (p_up .- p_do .+ rhog_MPa .* (z_up .- z_do)) ./ (p_up .- p_do);
    return flow_approx
end

function vc_integral(z_do, z_up, p_do, p_up, a, b, Kmax, flow_approx)
    denom = b *(Kmax * rhog_MPa *(z_up - z_do) * (a+1) / a .+ flow_approx);
    num_do = log.(Complex.(a * (Kmax * rhog_MPa *(z_up .- z_do) * (a+1) / a .+ flow_approx) .* exp.(b.*p_do) .+ flow_approx));
    num_up = log.(Complex.(a * (Kmax * rhog_MPa *(z_up .- z_do) * (a+1) / a .+ flow_approx) .* exp.(b.*p_up) .+ flow_approx));
    flow = -Kmax * (a+1)/ a * flow_approx .* (num_up .- num_do) ./ denom;
    return flow
end
 
function theta_to_p(theta)      
    p = (theta.-1).*5;  
    return p
end

function p_to_theta(p)
    theta = p./5 .+ 1; 
    return theta
end

function fstem!(F_stem,x)
    flow_in_stem_approx = vc_integral_approx(z_soil, z_bottom_stem, p_soil, x[1], a_root, b_root, K_max_root_total_moles)
    flow_in_stem = vc_integral(z_soil, z_bottom_stem, p_soil, x[1], a_root, b_root, K_max_root_total_moles, flow_in_stem_approx)  
    F_stem[1] = flow_in_stem - T_0   
end

function fleaf!(F_leaf,y)
    flow_out_stem_approx = vc_integral_approx(z_bottom_stem, z_leaf, p_stem_ini, y[1], a_stem, b_stem, K_max_stem_moles)
    flow_out_stem = vc_integral(z_bottom_stem, z_leaf, p_stem_ini, y[1], a_stem, b_stem, K_max_stem_moles, flow_out_stem_approx)    
    F_leaf[1] = flow_out_stem - T_0
end

# Solver
function roots(dy,y,paramset,t)
    p_stem = theta_to_p(y[1]/size_reservoir_stem_moles) 
    p_leaf = theta_to_p(y[2]/size_reservoir_leaf_moles) 

    flow_in_stem_approx = vc_integral_approx(z_soil, z_bottom_stem, p_soil, p_stem, a_root, b_root, K_max_root_total_moles)
    flow_in_stem = vc_integral(z_soil, z_bottom_stem, p_soil, p_stem, a_root, b_root, K_max_root_total_moles, flow_in_stem_approx)   
    flow_out_stem_approx = vc_integral_approx(z_bottom_stem, z_leaf, p_stem, p_leaf, a_stem, b_stem, K_max_stem_moles)
    flow_out_stem = vc_integral(z_bottom_stem, z_leaf, p_stem, p_leaf, a_stem, b_stem, K_max_stem_moles, flow_out_stem_approx) 
     
    if t < 500 
        T = T_0;
    elseif t >= 500 && t < 1000
        T = 10*(T_0/5)*(t-500)/500+T_0
    else t >= 1000
        T = 10*(T_0/5)*500/500+T_0
    end    
    dy[1] = flow_in_stem - flow_out_stem
    dy[2] = flow_out_stem - T
end

########## Initial values ###########
z_soil = 0.0
z_bottom_stem = 1.0
z_leaf = 14.0

#T_0 = 0.01/mass_mole_water # moles per s, using value at noon fig 9j) Christoffersen                    
#p_soil = -0.34 # MPa want it to be wet and wetter than rest of plant, so close to 0

T_0 = 0.0 
p_soil = -15
p_stem_ini = -0.5
p_leaf_ini = -0.05

#= Set system to equilibrium state by setting LHS of both odes to 0
solnstem = nlsolve(fstem!, [-1.0])
p_stem_ini = solnstem.zero[1]
solnleaf= nlsolve(fleaf!, [-1.0])
p_leaf_ini = solnleaf.zero[1]
=#

theta_stem_0 = p_to_theta(p_stem_ini) 
theta_leaf_0 = p_to_theta(p_leaf_ini) 
y1_0 = float(theta_stem_0*size_reservoir_stem_moles)
y2_0 = float(theta_leaf_0*size_reservoir_leaf_moles)
y0 = [y1_0; y2_0]

# Simulation length
tend = 60*60.0*2
tspan = (0.0,tend)
dt = 1
alg = Euler()
alg_name = "Euler"

# Solve the problem
prob = ODEProblem(roots,y0,tspan,paramset)
sol = solve(prob,alg,adaptive=false,dt=dt)

y_1 = reduce(hcat,sol.u)[1,:]
y_2 = reduce(hcat,sol.u)[2,:]

plot(sol.t, y_1, label="stem", xaxis="t [s]", yaxis="water content [mol]", dpi=500)
plot!(sol.t, y_2, label="leaf")
plot!(sol.t, y1_0.*ones(length(sol.t),1),label="W_stem_0",dpi=500)
plot!(sol.t, y2_0.*ones(length(sol.t),1),label="W_leaf_0",dpi=500)
savefig("water_content_moles.png") 

# Convert soln to volumetric water content and plot
y_theta_1 = y_1/size_reservoir_stem_moles
y_theta_2 = y_2/size_reservoir_leaf_moles

plot(sol.t, y_theta_1, label="stem", xaxis="t [s]", yaxis="relative water content [mol/mol]",dpi=500)
plot!(sol.t, y_theta_2, label="leaf", dpi=500)
plot!(sol.t, theta_stem_0.*ones(length(sol.t),1),label="w_stem_0",dpi=500)
plot!(sol.t, theta_leaf_0.*ones(length(sol.t),1),label="w_leaf_0",dpi=500)
savefig("relative_water_content.png") 

# Compute pressure and flow rates from soln and plot
p_stem = theta_to_p(y_theta_1) 
p_leaf = theta_to_p(y_theta_2) 

T = collect(0.0:dt:Float64(tend))
tlength = 1:length(T)
for i in tlength
if T[i] < 500
    T[i] = T_0;
elseif T[i] >= 500 && T[i] < 1000
    T[i] = 10*T_0/5*(T[i]-500)/500+T_0
else T[i] >= 1000
    T[i] = 10*(T_0/5)*500/500+T_0
end
end

# Plot pressure in stem and leaves as function of time [MPa]
plot(sol.t,p_stem,linewidth=2,xaxis="time [s]",yaxis="pressure [MPa]",label="stem",dpi=500)
plot!(sol.t,p_leaf,linewidth=2,label="leaf",dpi=500)
plot!(sol.t,p_stem_ini.*ones(length(sol.t),1),label="p_stem_0",dpi=500)
plot!(sol.t,p_leaf_ini.*ones(length(sol.t),1),label="p_leaf_0",dpi=500)
savefig("pressure.png") 

# Plot flow as a function of time [mol s-1]
flow_in_stem_approx = vc_integral_approx(z_soil, z_bottom_stem, p_soil, p_stem, a_root, b_root, K_max_root_total_moles)
flow_in_stem = vc_integral(z_soil, z_bottom_stem, p_soil, p_stem, a_root, b_root, K_max_root_total_moles, flow_in_stem_approx)   
flow_out_stem_approx = vc_integral_approx(z_bottom_stem, z_leaf, p_stem, p_leaf, a_stem, b_stem, K_max_stem_moles)
flow_out_stem = vc_integral(z_bottom_stem, z_leaf, p_stem, p_leaf, a_stem, b_stem, K_max_stem_moles, flow_out_stem_approx) 
plot(sol.t,real.(flow_in_stem),linewidth=2,xaxis="time [s]",yaxis="flow [mol s-1]",label="flow into stem",legend=:bottomright,dpi=500)
plot!(sol.t,real.(flow_out_stem),linewidth=2,label="flow into leaves",dpi=500)
plot!(sol.t,T,linewidth=2,label="transpiration boundary condition",dpi=500)

flow_in_stem_approx = vc_integral_approx(z_soil, z_bottom_stem, p_soil, p_stem_ini, a_root, b_root, K_max_root_total_moles)
flow_in_stem_ini = vc_integral(z_soil, z_bottom_stem, p_soil, p_stem_ini, a_root, b_root, K_max_root_total_moles, flow_in_stem_approx)  
flow_out_stem_approx = vc_integral_approx(z_bottom_stem, z_leaf, p_stem_ini, p_leaf_ini, a_stem, b_stem, K_max_stem_moles)
flow_out_stem_ini = vc_integral(z_bottom_stem, z_leaf, p_stem_ini, p_leaf_ini, a_stem, b_stem, K_max_stem_moles, flow_out_stem_approx)    
plot!(sol.t,real(flow_in_stem_ini).*ones(length(sol.t),1),linewidth=2,xaxis="time [s]",yaxis="flow [mol s-1]",label="initial flow into stem",legend=:bottomright,dpi=500)
plot!(sol.t,real(flow_out_stem_ini).*ones(length(sol.t),1),linewidth=2,label="initial flow into leaves",dpi=500)
savefig("flow.png") 
