using XLSX
using DifferentialEquations
using Plots
using NLsolve
using Revise

# Load parameter file
xf = XLSX.readxlsx("/Users/annagl/Documents/LandHydrology.jl/src/parameters_plant_hydraulics.xlsx")
XLSX.sheetnames(xf)
sh = xf["Sheet1"] 

# Parameters
rhog_MPa = sh[5,3]
mass_mole_water = sh[6,3]
h_root = sh[8,3]
h_stem = sh[9,3]
K_max_stem_moles = sh[43,3]
K_max_effective_root_moles = sh[45,3]
size_reservoir_stem_moles = sh[54,3]
size_reservoir_leaf_moles = sh[55,3]
a_root = sh[78,3]
a_stem = sh[79,3] 
b_root = sh[80,3]
b_stem = sh[81,3]

paramset = [rhog_MPa mass_mole_water h_root h_stem K_max_stem_moles K_max_effective_root_moles size_reservoir_stem_moles size_reservoir_leaf_moles a_root a_stem b_root b_stem];

# Functions
function vc_integral_approx(z_do, z_up, p_do, p_up, a, b, Kmax)
    u_do = a .* exp.(b .* p_do) 
    u_up = a .* exp.(b .* p_up) 
    num_do = log.(u_do .+ 1)
    num_up = log.(u_up .+ 1)
    c = Kmax .* (a+1) / a
    d = rhog_MPa .* (z_up .- z_do)
    flow_approx = -c ./ b .* (num_up .- num_do) .* (p_up .- p_do .+ d) ./ (p_up .- p_do)
    A = c .* d .+ flow_approx
    B = -c .* flow_approx ./ (b .* A)
    return u_do, u_up, A, B, flow_approx
end

function vc_integral(u_do, u_up, A, B, flow_approx)   
    flow = B .* (log.(u_up .* A .+ flow_approx) .- log.(u_do .* A .+ flow_approx))
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
    u_do, u_up, A, B, flow_in_stem_approx = vc_integral_approx(z_soil, z_bottom_stem, p_soil, x[1], a_root, b_root, K_max_effective_root_moles)
    flow_in_stem = vc_integral(u_do, u_up, A, B, flow_in_stem_approx)  
    F_stem[1] = flow_in_stem - T_0   
end

function fleaf!(F_leaf,y)
    u_do, u_up, A, B, flow_out_stem_approx = vc_integral_approx(z_bottom_stem, z_leaf, p_stem_ini, y[1], a_stem, b_stem, K_max_stem_moles)
    flow_out_stem = vc_integral(u_do, u_up, A, B, flow_out_stem_approx) 
    F_leaf[1] = flow_out_stem - T_0
end

# Solver
function roots(dy,y,paramset,t)
    p_stem = theta_to_p(y[1]/size_reservoir_stem_moles) 
    p_leaf = theta_to_p(y[2]/size_reservoir_leaf_moles) 

    u_do, u_up, A, B, flow_in_stem_approx = vc_integral_approx(z_soil, z_bottom_stem, p_soil, p_stem, a_root, b_root, K_max_effective_root_moles)
    flow_in_stem = vc_integral(u_do, u_up, A, B, flow_in_stem_approx) 
    u_do, u_up, A, B, flow_out_stem_approx = vc_integral_approx(z_bottom_stem, z_leaf, p_stem, p_leaf, a_stem, b_stem, K_max_stem_moles)
    flow_out_stem = vc_integral(u_do, u_up, A, B, flow_out_stem_approx) 
    
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

# Initial values
T_0 = 0.01/mass_mole_water # moles per s, using value at noon fig 9j) Christoffersen                    
p_soil = -0.02 # MPa
z_soil = 0.0 # m
z_bottom_stem = h_root # m
z_leaf = h_stem # m

# Set system to equilibrium state by setting LHS of both ODEs to 0
solnstem = nlsolve(fstem!, [-1.0])
p_stem_ini = solnstem.zero[1]
solnleaf= nlsolve(fleaf!, [-1.0])
p_leaf_ini = solnleaf.zero[1]

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
savefig("absolute_water_content.png") 

# Convert soln to volumetric water content and plot
y_theta_1 = y_1/size_reservoir_stem_moles
y_theta_2 = y_2/size_reservoir_leaf_moles
plot(sol.t, y_theta_1, label="stem", xaxis="t [s]", yaxis="relative water content [m3/m3]",dpi=500)
plot!(sol.t, y_theta_2, label="leaf", dpi=500)
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
u_do, u_up, A, B, flow_in_stem_approx = vc_integral_approx(z_soil, z_bottom_stem, p_soil, p_stem, a_root, b_root, K_max_effective_root_moles)
flow_in_stem = vc_integral(u_do, u_up, A, B, flow_in_stem_approx) 
u_do, u_up, A, B, flow_out_stem_approx = vc_integral_approx(z_bottom_stem, z_leaf, p_stem, p_leaf, a_stem, b_stem, K_max_stem_moles)
flow_out_stem = vc_integral(u_do, u_up, A, B, flow_out_stem_approx) 
plot(sol.t,flow_in_stem,linewidth=2,xaxis="time [s]",yaxis="flow [mol s-1]",label="flow into stem",legend=:bottomright,dpi=500)
plot!(sol.t,flow_out_stem,linewidth=2,label="flow into leaves",dpi=500)
plot!(sol.t,T,linewidth=2,label="transpiration boundary condition",dpi=500)
#savefig("flow.png") 