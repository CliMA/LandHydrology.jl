using XLSX
using DifferentialEquations
using Plots
using NLsolve
using Revise

# Load parameter file
xf = XLSX.readxlsx("/Users/annagl/Documents/LandHydrology.jl/src/parameters_plant_hydraulics_original.xlsx")
XLSX.sheetnames(xf)
sh = xf["Sheet1"] 

# Parameters
rhog_MPa = sh[5,3]
mass_mole_water = sh[6,3]
h_root = sh[8,3]
h_stem = sh[9,3]
K_max_stem_moles = sh[43,3]*h_stem
K_max_effective_root_moles = sh[45,3]*h_root
K_max_leaf_moles = K_max_stem_moles*1.3
size_reservoir_stem_moles = sh[54,3]
size_reservoir_leaf_moles = sh[55,3]
a_root = sh[78,3]
a_stem = sh[79,3] 
a_leaf = a_stem*1.3
b_root = sh[80,3]
b_stem = sh[81,3]
b_leaf = b_stem*1.3

pi_0 = sh[63,3]
pi_tlp = sh[65,3]
theta_r = sh[67,3]
theta_tlp = sh[69,3]
theta_ft = sh[72,3]
epsilon = sh[61,3]
S_s = sh[60,3]

paramset = [rhog_MPa mass_mole_water h_root h_stem K_max_stem_moles K_max_effective_root_moles K_max_leaf_moles size_reservoir_stem_moles size_reservoir_leaf_moles a_root a_stem a_leaf b_root b_stem b_leaf pi_0 pi_tlp theta_r theta_tlp theta_ft epsilon S_s] 

# Functions
function flow(z_do, z_up, p_do, p_up, a, b, Kmax, a_up, b_up, Kmax_up)
    u_do = a .* exp.(b .* p_do) 
    u_up = a .* exp.(b .* p_up) 
    num_do = log.(u_do .+ 1)
    num_up = log.(u_up .+ 1)
    c = Kmax .* (a+1) / a
    term1 = -c ./ b .* (num_up .- num_do) ./(z_up .- z_do)

    c_up = Kmax_up .* (a_up .+ 1) ./ a_up
    term2_up = -c_up .* (1 .- 1 ./ (1 .+ a_up.*exp.(b_up .* p_up))) .* rhog_MPa
    term2_do = -c .* (1 .- 1 ./ (1 .+ a.*exp.(b .* p_do))) .* rhog_MPa
    term2 = (term2_up .+ term2_do)./2
    flow = term1 .+ term2
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

#=

function theta_to_p(theta, pi_0, theta_r, theta_tlp, theta_ft, epsilon, S_s) 
    if theta > theta_ft      
        p = ((theta-theta_ft)./S_s)^2;
    elseif theta > theta_tlp && theta <= theta_ft
        p_sol = -abs(pi_0).*(theta_ft-theta_r)./(theta-theta_r);
        p_p = abs(pi_0)-epsilon.*(theta_ft-theta)./(theta_ft-theta_r);
        p = p_sol + p_p; 

    elseif theta >= theta_r && theta <= theta_tlp
        p_sol = -abs(pi_0).*(theta_ft-theta_r)./(theta-theta_r);
        p = p_sol;
    else
        @show("theta out of bounds")
    end 
    return p  
end

function p_to_theta(p, pi_0, pi_tlp, theta_r, theta_ft, epsilon, S_s)
    if p > 0
        theta = sqrt(p).*S_s+theta_ft;
    elseif p > pi_tlp && p <= 0
        theta = (epsilon*theta_r + epsilon*theta_ft - p.*theta_r + p.*theta_ft + theta_r*abs(pi_0) - theta_ft*abs(pi_0) - theta_r*(abs(pi_0).^2 + 2*epsilon.*p + 2*epsilon*abs(pi_0) - 2*p.*abs(pi_0) + epsilon^2 + p.^2).^(1/2) + theta_ft*(abs(pi_0)^2 + 2*epsilon.*p + 2*epsilon*abs(pi_0) - 2*p*abs(pi_0) + epsilon^2 + p.^2).^(1/2))/(2*epsilon);
    elseif p <= pi_tlp
        theta = theta_r + (abs(pi_0)*(theta_r - theta_ft))./p;
    else 
    @show("p out of bounds")
    end 
    return theta  
end
=#

function fstem!(F_stem,x)
    flow_in_stem = flow(z_soil, z_bottom_stem, p_soil, x[1], a_root, b_root, K_max_effective_root_moles, a_stem, b_stem, K_max_stem_moles) 
    F_stem[1] = flow_in_stem - T_0   
end

function fleaf!(F_leaf,y)
    flow_out_stem = flow(z_bottom_stem, z_leaf, p_stem_ini, y[1], a_stem, b_stem, K_max_stem_moles, a_leaf, b_leaf, K_max_leaf_moles)
    F_leaf[1] = flow_out_stem - T_0
end

# Solver
function roots(dy,y,paramset,t)
    p_stem = theta_to_p(y[1]/size_reservoir_stem_moles) 
    p_leaf = theta_to_p(y[2]/size_reservoir_leaf_moles) 

    #p_stem = theta_to_p(y[1]/size_reservoir_stem_moles, pi_0, theta_r, theta_tlp, theta_ft, epsilon, S_s)
    #p_leaf = theta_to_p(y[2]/size_reservoir_leaf_moles, pi_0, theta_r, theta_tlp, theta_ft, epsilon, S_s)

    flow_in_stem = ones(n_roots,1)
    for i in 1:n_roots
        flow_in_stem[i] = flow(z_soil[i], z_bottom_stem, p_soil[i], p_stem, a_root, b_root, K_max_effective_root_moles, a_stem, b_stem, K_max_stem_moles)
        #@show(flow_in_stem)
    end

    #flow_in_stem = flow(z_soil, z_bottom_stem, p_soil, p_stem, a_root, b_root, K_max_effective_root_moles, a_stem, b_stem, K_max_stem_moles)
    flow_out_stem = flow(z_bottom_stem, z_leaf, p_stem, p_leaf, a_stem, b_stem, K_max_stem_moles, a_leaf, b_leaf, K_max_leaf_moles)
    
    if t < 500 
        T = T_0;
    elseif t >= 500 && t < 1000
        T = 10*(T_0/5)*(t-500)/500+T_0
    else t >= 1000
        T = 10*(T_0/5)*500/500+T_0
    end

    dy[1] = sum(flow_in_stem) - flow_out_stem
    dy[2] = flow_out_stem - T
end

# Initial values
T_0 = 0.01/mass_mole_water # moles per s, using value at noon fig 9j) Christoffersen                    
p_soil = [-0.1; -0.2; -0.3] # MPa want it to be wet and wetter than rest of plant, so close to 0
z_soil = [0; 1; 2]
n_roots = length(p_soil)

#p_soil = -0.02 # MPa
#z_soil = 0.0 # m
z_bottom_stem = 5 #h_root # m
z_leaf = h_stem # m

# Set system to equilibrium state by setting LHS of both ODEs to 0
#solnstem = nlsolve(fstem!, [-1.0])
p_stem_ini = -0.5 #solnstem.zero[1]
#solnleaf= nlsolve(fleaf!, [-1.0])
p_leaf_ini = -0.4 #solnleaf.zero[1]

theta_stem_0 = p_to_theta(p_stem_ini) 
theta_leaf_0 = p_to_theta(p_leaf_ini) 

#theta_stem_0 = p_to_theta(p_stem_ini, pi_0, pi_tlp, theta_r, theta_ft, epsilon, S_s)
#theta_leaf_0 = p_to_theta(p_leaf_ini, pi_0, pi_tlp, theta_r, theta_ft, epsilon, S_s)

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

plot(sol.t, y_1, label="stem", xaxis="t [s]", yaxis="absolute water content [mol]", dpi=500)
plot!(sol.t, y_2, label="leaf")
plot!(sol.t, y1_0.*ones(length(sol.t),1),label="W_stem_0",dpi=500)
plot!(sol.t, y2_0.*ones(length(sol.t),1),label="W_leaf_0",dpi=500)
savefig("absolute_water_content2.0.png") 

# Convert soln to volumetric water content and plot
y_theta_1 = y_1/size_reservoir_stem_moles
y_theta_2 = y_2/size_reservoir_leaf_moles
plot(sol.t, y_theta_1, label="stem", xaxis="t [s]", yaxis="relative water content [m3/m3]",dpi=500)
plot!(sol.t, y_theta_2, label="leaf", dpi=500)
savefig("relative_water_content2.0.png") 

# Compute pressure and flow rates from soln and plot
p_stem = theta_to_p(y_theta_1) 
p_leaf = theta_to_p(y_theta_2) 

#p_stem = theta_to_p.(y_theta_1, pi_0, theta_r, theta_tlp, theta_ft, epsilon, S_s)
#p_leaf = theta_to_p.(y_theta_2, pi_0, theta_r, theta_tlp, theta_ft, epsilon, S_s)

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
savefig("pressure2.0.png") 

# Plot flow as a function of time [mol s-1]
flow_in_stem = ones(length(p_stem),n_roots)
for j in 1:length(p_stem)
    for i in 1:n_roots
        flow_in_stem[j,i] = flow(z_soil[i], z_bottom_stem, p_soil[i], p_stem[j], a_root, b_root, K_max_effective_root_moles, a_stem, b_stem, K_max_stem_moles)
    end
end

#flow_in_stem = flow(z_soil, z_bottom_stem, p_soil, p_stem, a_root, b_root, K_max_effective_root_moles, a_stem, b_stem, K_max_stem_moles)
flow_in_stem_root_1 = flow_in_stem[:,1]
flow_in_stem_root_2 = flow_in_stem[:,2]
flow_in_stem_root_3 = flow_in_stem[:,3]
total_flow_in_stem = sum(flow_in_stem, dims=2)
flow_out_stem = flow(z_bottom_stem, z_leaf, p_stem, p_leaf, a_stem, b_stem, K_max_stem_moles, a_leaf, b_leaf, K_max_leaf_moles)
plot(sol.t,flow_in_stem_root_1,linewidth=2,xaxis="time [s]",yaxis="flow [mol s-1]",label="flow into stem from root 1",legend=:bottomright,dpi=500)
plot!(sol.t,flow_in_stem_root_2,linewidth=2,xaxis="time [s]",yaxis="flow [mol s-1]",label="flow into stem from root 2",legend=:bottomright,dpi=500)
plot!(sol.t,flow_in_stem_root_3,linewidth=2,xaxis="time [s]",yaxis="flow [mol s-1]",label="flow into stem from root 3",legend=:bottomright,dpi=500)
plot!(sol.t,total_flow_in_stem,linewidth=2,xaxis="time [s]",yaxis="flow [mol s-1]",label="total flow into stem",legend=:bottomright,dpi=500)
plot!(sol.t,flow_out_stem,linewidth=2,label="flow into leaves",dpi=500)
plot!(sol.t,T,linewidth=2,label="transpiration boundary condition",dpi=500)
savefig("flow2.0.png") 