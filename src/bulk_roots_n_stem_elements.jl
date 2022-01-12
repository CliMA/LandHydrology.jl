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

# Set number of stem elements
n_stem_elements = 2

# Set number of leaves
n_leaves = 2
n_tot = n_stem_elements + n_leaves 
paramset = [n_stem_elements n_leaves a_root a_stem b_root b_stem K_max_root_total_moles K_max_stem_moles rho_water rhog_MPa h_root h_stem pi_0_leaf pi_0_stem epsilon_leaf epsilon_stem theta_r_leaf theta_r_stem WD theta_tlp_leaf theta_tlp_stem theta_ft_stem size_reservoir_stem_moles size_reservoir_leaf_moles theta_ft_leaf volume_mole_water S_stem S_leaf porosity_stem porosity_leaf];

########## Functions #########
 
    function vc_integral_approx(p_dos, p_ups, a, b)
    #     if p_ups > 0 && p_dos > 0
    #         error('positive pressure value in approx integral')
    #     end
        t_dos = log.(a .* exp.(b.*p_dos) .+ 1);
        t_ups = log.(a .* exp.(b.*p_ups) .+ 1);
        flux_approx = (a+1) / (a*b) .* (t_ups .- t_dos);
        return flux_approx
    end
    
    function vc_integral(p_dos, p_ups, h, flux_approx, Kmax, rhog_MPa, a, b)
    #     if p_ups <= 0 && p_dos <= 0
    #         error('positive pressure value in approx integral')
    #     end
        krghe = Kmax * rhog_MPa * h * (a+1) / a .+ flux_approx;
        lower = b .* krghe;
        multi = Kmax * (a+1) / a .* flux_approx;
        upper_dos = log.(a .* krghe .* exp.(b.*p_dos) .+ flux_approx);
        upper_ups = log.(a .* krghe .* exp.(b.*p_ups) .+ flux_approx);
        flux = multi .* (upper_ups .- upper_dos) ./ lower;
        return flux
    end
     
    function theta_to_p(theta, pi_0, theta_r, theta_tlp, theta_ft, epsilon, S_s)      
        p = (theta.-1).*5;  
        return p
    end
    
    function p_to_theta(p, pi_0, pi_tlp, theta_r, theta_ft, epsilon, S_s)
        theta = p./5 .+ 1; 
        return theta
    end
    
    #= function p = theta_to_p(theta, pi_0, theta_r, theta_tlp, theta_ft, epsilon, S_s) 
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
             error('theta out of bounds')
         end   
    end
    function [theta] = p_to_theta(p, pi_0, pi_tlp, theta_r, theta_ft, epsilon, S_s)
         if p > 0
             theta = sqrt(p).*S_s+theta_ft;
         elseif p > pi_tlp && p <= 0
             theta = (epsilon*theta_r + epsilon*theta_ft - p.*theta_r + p.*theta_ft + theta_r*abs(pi_0) - theta_ft*abs(pi_0) - theta_r*(abs(pi_0).^2 + 2*epsilon.*p + 2*epsilon*abs(pi_0) - 2.*p.*abs(pi_0) + epsilon^2 + p.^2).^(1/2) + theta_ft*(abs(pi_0)^2 + 2*epsilon.*p + 2*epsilon*abs(pi_0) - 2.*p*abs(pi_0) + epsilon^2 + p.^2).^(1/2))/(2*epsilon);
         elseif p <= pi_tlp
             theta = theta_r + (abs(pi_0)*(theta_r - theta_ft))./p;
         else 
             error('p out of bounds')
         end   
    end
    =#

function fstem!(F_stem,x)   
    flow = ones(1,n_stem_elements)
    flow_in_stem_approx = K_max_root_total_moles*vc_integral_approx(x[1], p_soil, a_root, b_root)*(p_soil - x[1] - rhog_MPa* h_root/2) / (p_soil - x[1])
    flow[1] = vc_integral(x[1], p_soil, h_root/2, flow_in_stem_approx,K_max_root_total_moles, rhog_MPa, a_root, b_root)
    
    if n_stem_elements > 1
        for i in 2:n_stem_elements-1
            flow_in_stem_approx = K_max_stem_moles*vc_integral_approx(x[i], x[i-1], a_root, b_root)*(x[i-1] - x[i] - rhog_MPa* h_stem/n_stem_elements/2) / (x[i-1] - x[i])
            flow[i] = vc_integral(x[i], x[i-1], h_stem/n_stem_elements/2, flow_in_stem_approx, K_max_stem_moles, rhog_MPa, a_root, b_root)   
            F_stem[i] = flow[i-1] - flow[i]
        end
    end

    F_stem[end] = flow[end] - T_0   
end

function fleaf!(F_leaf,y)
    # Leaf 1
    flow_in_leaf_approx = K_max_stem_moles*vc_integral_approx(y[1], p_stem_ini[end], a_stem, b_stem)*(p_stem_ini[end] - y[1] - rhog_MPa* h_stem/n_stem_elements/2) / (p_stem_ini[end] - y[1])
    flow[n_stem_elements+1] = vc_integral(y[1], p_stem[end], h_stem/n_stem_elements/2, flow_in_leaf_approx, K_max_stem_moles, rhog_MPa, a_stem, b_stem)  

    # Leaf n
    for i in 2:n_leaves
        flow_in_leaf_approx = K_max_stem_moles*vc_integral_approx(p_leaf[i], p_leaf[i-1], a_stem, b_stem)*(p_leaf[i-1] - p_leaf[i] - rhog_MPa* h_stem/n_stem_elements/2) / (p_leaf[i-1] - p_leaf[i])
        flow[n_stem_elements+2] = vc_integral(p_leaf[i], p_leaf[i-1], h_stem/n_stem_elements/2, flow_in_stem_approx, K_max_stem_moles, rhog_MPa, a_stem, b_stem)  
        @show(flow)
    end
end

########## Solver ###########

function roots(dy,y,paramset,t)
    
    p_stem = theta_to_p(y[1:n_stem_elements]./size_reservoir_stem_moles, pi_0_stem, theta_r_stem, theta_tlp_stem, theta_ft_stem, epsilon_stem, S_stem)      
    p_leaf = theta_to_p(y[n_stem_elements+1:end]./size_reservoir_leaf_moles, pi_0_leaf, theta_r_leaf, theta_tlp_leaf, theta_ft_leaf, epsilon_leaf, S_leaf)
    flow = ones(1,n_stem_elements+n_leaves)
    @show(p_stem)
    @show(p_leaf)

    # Stem 1
    flow_in_stem_approx = K_max_root_total_moles*vc_integral_approx(p_stem[1], p_soil, a_root, b_root)*(p_soil - p_stem[1] - rhog_MPa * h_root/2) / (p_soil - p_stem[1])
    flow[1] = vc_integral(p_stem[1], p_soil, h_root/2, flow_in_stem_approx, K_max_root_total_moles, rhog_MPa, a_root, b_root)   

    # Stem n
    for i in 2:n_stem_elements
        @show(i)
        flow_in_stem_approx = K_max_stem_moles*vc_integral_approx(p_stem[i], p_stem[i-1], a_root, b_root)*(p_stem[i-1] - p_stem[i] - rhog_MPa* h_stem/n_stem_elements/2) / (p_stem[i-1] - p_stem[i])
        flow[i] = vc_integral(p_stem[i], p_stem[i-1], h_stem/n_stem_elements/2, flow_in_stem_approx,K_max_stem_moles, rhog_MPa, a_root, b_root)   
        @show(flow)
    end

    # Leaf 1
    flow_in_leaf_approx = K_max_stem_moles*vc_integral_approx(p_leaf[1], p_stem[end], a_stem, b_stem)*(p_stem[end] - p_leaf[1] - rhog_MPa* h_stem/n_stem_elements/2) / (p_stem[end] - p_leaf[1])
    flow[n_stem_elements+1] = vc_integral(p_leaf[1], p_stem[end], h_stem/n_stem_elements/2, flow_in_leaf_approx, K_max_stem_moles, rhog_MPa, a_stem, b_stem)  

    # Leaf n
    for i in 2:n_leaves
        flow_in_leaf_approx = K_max_stem_moles*vc_integral_approx(p_leaf[i], p_leaf[i-1], a_stem, b_stem)*(p_leaf[i-1] - p_leaf[i] - rhog_MPa* h_stem/n_stem_elements/2) / (p_leaf[i-1] - p_leaf[i])
        flow[n_stem_elements+2] = vc_integral(p_leaf[i], p_leaf[i-1], h_stem/n_stem_elements/2, flow_in_stem_approx, K_max_stem_moles, rhog_MPa, a_stem, b_stem)  
        @show(flow)
    end

    if t < 500 
        T = T_0;
    elseif t >= 500 && t < 1000
        T = 10*(T_0/5)*(t-500)/500+T_0
    else t >= 1000
        T = 10*(T_0/5)*500/500+T_0
    end
    
    for i in 1:n_tot-1
        dy[i] = flow[i] - flow[i+1]
    end

    dy[end] = flow[end] - T
end

########## Initial values ###########
T_0 = 0.01/mass_mole_water # moles per s, using value at noon fig 9j) Christoffersen                    
p_soil = 0.0 # MPa want it to be saturated, and more saturated than rest of plant, so close to 0

# Set system to equilibrium state by setting LHS of both odes to 0
starting_point = -1.0

solnstem = nlsolve(fstem!, starting_point.*ones(n_stem_elements,1))
p_stem_1_ini = solnstem.zero[1]

#=
solnleaf= nlsolve(fleaf!, [starting_point])
p_leaf_2_ini = solnleaf.zero[1]

theta_stem_0 = p_to_theta(p_stem_ini, pi_0_stem, pi_tlp_stem, theta_r_stem, theta_ft_stem, epsilon_stem, S_stem)
theta_leaf_0 = p_to_theta(p_leaf_ini, pi_0_leaf, pi_tlp_leaf, theta_r_leaf, theta_ft_leaf, epsilon_leaf, S_leaf)

y1_0 = float(theta_stem_0*size_reservoir_stem_moles) #.*(ones(n_stem_elements,1))
y2_0 = y1_0*1.2
y3_0 = float(theta_leaf_0*size_reservoir_leaf_moles) #.*(ones(n_leaves,1))
y4_0 = y3_0*1.2

y0 = [y1_0; y2_0; y3_0; y4_0]

############### Simulation length ###
tend = 60*60.0*2
tspan = (0.0,tend)
dt = 0.1
alg = Euler()
alg_name = "Euler"

#alg = GIRK4()
#alg = SSPRK22()
#alg_name = "SSPRK22"

################ Solve the problem

prob = ODEProblem(roots,y0,tspan,paramset)
sol = solve(prob,alg,adaptive=false,dt=dt)

for i in 1:n_stem_elements
    plot(sol.t, reduce(hcat,sol.u)[i,:], label="stem", xaxis="t [s]", yaxis="water content [mol]", dpi=500)
end 

for i in n_stem_elements+1:n_tot
    plot!(sol.t, reduce(hcat,sol.u)[i,:], label="leaves")
end 
savefig("water_content_moles.png")

# Convert soln to volumetric water content and plot
y_theta_1 = y_1/size_reservoir_stem_moles
y_theta_2 = y_2/size_reservoir_leaf_moles

plot(sol.t,y_theta_1, label="stem", xaxis="t [s]", yaxis="relative water content [mol/mol]",dpi=500)
plot!(sol.t,y_theta_2, label="leaves",dpi=500)
savefig("relative_water_content.png") 

# Compute pressure and flow rates from soln and plot
p_stem = theta_to_p(y_theta_1, pi_0_stem, theta_r_stem, theta_tlp_stem, theta_ft_stem, epsilon_stem, S_stem)       
p_leaf = theta_to_p(y_theta_2, pi_0_leaf, theta_r_leaf, theta_tlp_leaf, theta_ft_leaf, epsilon_leaf, S_leaf)
flow_in_stem_approx = K_max_root_total_moles.*vc_integral_approx(p_stem, p_soil, a_root, b_root).*(p_soil .- p_stem .- rhog_MPa * h_root/2) ./ (p_soil .- p_stem)
flow_in_stem = vc_integral(p_stem, p_soil, h_root/2, flow_in_stem_approx,K_max_root_total_moles, rhog_MPa, a_root, b_root)   
flow_out_stem_approx = K_max_stem_moles.*vc_integral_approx(p_leaf, p_stem, a_stem, b_stem).*(p_stem .- p_leaf .- rhog_MPa*h_stem/2) ./ (p_stem .- p_leaf)
flow_out_stem = vc_integral(p_leaf, p_stem, h_stem/2, flow_out_stem_approx, K_max_stem_moles, rhog_MPa, a_stem, b_stem)  

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
plot!(sol.t,p_leaf,linewidth=2,label="leaves",dpi=500)
savefig("pressure.png") 

# Plot flow as a function of time [mol s-1]
plot(sol.t,flow_in_stem,linewidth=2,xaxis="time [s]",yaxis="flow [mol s-1]",label="flow into stem",legend=:bottomright,dpi=500)
plot!(sol.t,flow_out_stem,linewidth=2,label="flow into leaves",dpi=500)
plot!(sol.t,T,linewidth=2,label="transpiration boundary condition",dpi=500)

savefig("flow.png") 

# Does this verify water conservation or just the implementation of the method? We are comparing change in 
# total water content over one time step to the net flow... I would say yes
# this is a way of verifying water cons as we are doing flow in stem - T (not looking at flow out stem)

# Plot error as function of time [mol]
yt = y_1 .+ y_2
ychange = diff(yt)
su = (flow_in_stem.-T).*dt
water_cons = ychange.-su[1:end-1] # 1:end-1 not 2:end because in Euler method y(t+h)=y(t0)+h*y'(t0), and t0 is time step 1
tvec = collect(1:length(ychange)) # We expect this to be very close to 0, because we are here implementing the Euler 1 method, could be error due to subtraction?
# interesting, when we do : water_cons = ychange.-su[1:end-1], we get exactly twice the error due to Euler method (2* 0.555 = 0.111 moles) 
plot(tvec,water_cons,xaxis="time [s]",yaxis="water conservation error per step [mol]", label=alg_name, title="diff(yt)-[(flow in stem-T).*dt](1:end-1)", dpi=500)
savefig("water_conservation_per_step.png") 

# Cumulative error over time
yt0 = y_1[1] .+ y_2[1]
water_cons_cum = (yt.-yt0.-(cumsum(su)))
plot(sol.t,water_cons_cum,xaxis="t",yaxis="cumulative water conservation error [mol]", label=alg_name, title="yt-y0-(cumsum((flow in stem-T).*dt))", dpi=500)
savefig("water_conservation_cumulative.png") 

# Explain error Forward Euler (Euler1)
sol_FE_1 = yt
eval_points_t2 = collect(0.0:dt:Float64(tend))
eval_points_t1 = zeros(length(eval_points_t2),1)
sol_exact = sum(y0) .+ analytical_approx(T_0,eval_points_t1,eval_points_t2)
plot(sol.t, sol_FE_1, label="FE soln",dpi=500)
plot!(sol.t, sol_exact, label="approximated analytic soln",dpi=500)
savefig("numerical_vs_approximated_analytic_soln.png") 

# Check water conservation with analyitc soln
plot(sol.t, sol_FE_1-sol_exact,dpi=500)
#savefig("error.png") 

# Compare to analytic soln
#=
l2norm = sqrt(sum((matlab_soln .- sol.u)^2))
plot(sol,linewidth=2,xaxis="t",label=["θ [rad]" "ω [rad/s]"],layout=(2,1))
=#

# Error : compare to approx

# Setting T=0, initial pressures in pools to 0 (saturation), and pressure in roots to something very negative (very dry) and see if the pools drain to soil, so negative flow rates. Pools should eventually empty out completely

# Setting T=0, pressure in roots to 0, pools to very negative, and see if they fill up, should see positive flow rates

# Plot with increasing number of compartments
# aim to show convergence

# Plot with increasing number of roots, different depths

# Compare model to data

# Compare model to plot in a paper

# Hook up to soil

# Soil model, lab experiment (?)

# Ask Yujie about tests he's done with data

# Change the boundary condition to a pressure one (VPD) 

# have a default grass, default plant, default tree ; different methods for each type once we start generalizing
# use similar structure as for soil, think about hierarchy and regrouping

P = linspace(-15,0,1000);
k_stem_power = (1 + (P./P_50_stem).^(a_x_stem)).^-1;
logisticEqn = '(a_stem+1)/a_stem*(1-1/(1+a_stem*exp(b_stem*x)))';
startPoints = [1 1];
fstem = fit(P',k_stem_power', logisticEqn,'Start', startPoints); 
stem_coeffs = coeffvalues(fstem);
a_stem = stem_coeffs(1);
b_stem = stem_coeffs(2);

k_root_power = (1 + (P./P_50_root).^(a_x_root)).^-1;
logisticEqn = '(a_root+1)/a_root*(1-1/(1+a_root*exp(b_root*x)))';
startPoints = [1 1];
froot = fit(P',k_root_power', logisticEqn,'Start', startPoints); 
root_coeffs = coeffvalues(froot);
a_root = root_coeffs(1);
b_root = root_coeffs(2);

=#