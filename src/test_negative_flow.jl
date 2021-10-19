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

    # p_dos must be dryier (more negative, so smaller) than p_ups
    
    #=
    p_dos = collect(-1:0.1:-0.5)
    p_ups = collect(-0.5:0.1:0.0)
    @show(p_dos .< p_ups) # p_dos is always smaller
    t_dos = log.(2 .* exp.(10 .*p_dos) .+ 1);
    t_ups = log.(2 .* exp.(10 .*p_ups) .+ 1);
    @show(t_dos .< t_ups)
    flux_approx = (2+1) / (2*10) .* (t_ups .- t_dos); # the soln is always positive

    plot(p_ups, t_ups)
    plot!(p_dos, t_dos)
    plot(collect(1:1:length(flux_approx)), flux_approx)
    =#

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

    # The term upriver (or upstream) refers to the 
    # direction towards the source of the river, i.e. against the direction of flow. Likewise, the term downriver (or downstream) describes the direction towards the mouth of the river, in which the current flows. 

    # We can have different scenarios :

    # Int w(p)*dP = P_ups - P_dos - rho * g * h
    
    # Pressure in the roots is more negative (smaller, dryer) than pressure in the stem
    # Water should flow from stem to roots
    # Sign of flow should be negative
    # Sign of h is positive
    # Sign of g is ...

    # Pressure in the stem is more negative (smaller, dryer) than pressure in the roots
    # Water should flow from stem to roots
    # Sign of flow should be positive
    # Sign of h is negative
    # Sign of g is ...

    function theta_to_p(theta, pi_0, theta_r, theta_tlp, theta_ft, epsilon, S_s)      
        p = (theta.-1).*5;  
        return p
    end
    
    function p_to_theta(p, pi_0, pi_tlp, theta_r, theta_ft, epsilon, S_s)
        theta = p./5 .+ 1; 
        return theta
    end

function fbase!(F_base,x)
    flow_in_base_approx = K_max_root_total_moles*vc_integral_approx(x[1], p_soil, a_root, b_root)*(p_soil - x[1] - rhog_MPa* h_root/2) / (p_soil - x[1])
    flow_in_base = vc_integral(x[1], p_soil, h_root/2, flow_in_base_approx,K_max_root_total_moles, rhog_MPa, a_root, b_root)
    F_base[1] = flow_in_base - T_0   
end

function fleaf!(F_leaf,y)
    flow_out_base_approx = K_max_stem_moles*vc_integral_approx(y[1], p_base_ini, a_stem, b_stem)*(p_base_ini - y[1] - rhog_MPa* h_stem/2) / (p_base_ini - y[1])
    flow_out_base = vc_integral(y[1], p_base_ini, h_stem/2, flow_out_base_approx, K_max_stem_moles, rhog_MPa, a_stem, b_stem)
    F_leaf[1] = flow_out_base - T_0
end

########## Solver ###########
function roots(dy,y,paramset,t)
    p_base = theta_to_p(y[1]/size_reservoir_stem_moles, pi_0_stem, theta_r_stem, theta_tlp_stem, theta_ft_stem, epsilon_stem, S_stem)      
    p_leaf = theta_to_p(y[2]/size_reservoir_leaf_moles, pi_0_leaf, theta_r_leaf, theta_tlp_leaf, theta_ft_leaf, epsilon_leaf, S_leaf)
    
    if p_base < p_soil 
        p_dos = p_base
        p_ups = p_soil
        sign_p = 1
        pdif = sign_p*(p_ups - p_dos)
        grav = -1*rhog_MPa* h_root/2
        sumf = pdif + grav
        if sumf >= 0
            sign_flow = 1
        else
            sign_flow = -1
        end
    else
        p_dos = p_soil
        p_ups = p_base
        sign_p = -1
        pdif = sign_p*(p_ups - p_dos)
        grav = -1*rhog_MPa* h_root/2
        sumf = pdif + grav
        sign_flow = -1
    end
    
    flow_in_base_approx = K_max_root_total_moles*vc_integral_approx(p_dos, p_ups, a_root, b_root)*(p_ups - p_dos - rhog_MPa* sign_p* h_root/2) / (p_ups - p_dos)
    flow_in_base = sign_flow*vc_integral(p_dos, p_ups, sign_p*h_root/2, flow_in_base_approx, K_max_root_total_moles, rhog_MPa, a_root, b_root)   

    if p_leaf < p_base
        p_dos = p_leaf
        p_ups = p_base
        sign_p = 1
        pdif = sign_p*(p_ups - p_dos)
        grav = -1*rhog_MPa* h_stem/2
        sumf = pdif + grav
        if sumf >= 0
            sign_flow = 1
        else
            sign_flow = -1
        end
    else
        p_dos = p_base
        p_ups = p_leaf
        sign_p = -1
        pdif = sign_p*(p_ups - p_dos)
        grav = -1*rhog_MPa* h_stem/2
        sumf = pdif + grav
        sign_flow = -1
    end

    flow_out_base_approx = K_max_stem_moles*vc_integral_approx(p_dos, p_ups, a_stem, b_stem)*(p_ups - p_dos - rhog_MPa* sign_p*h_stem/2) / (p_ups - p_dos)
    flow_out_base = sign_flow*vc_integral(p_dos, p_ups, sign_p*h_stem/2, flow_out_base_approx, K_max_stem_moles, rhog_MPa, a_stem, b_stem)  

    T = T_0
    
    dy[1] = flow_in_base - flow_out_base
    dy[2] = flow_out_base - T
end

########## Initial values ###########
T_0 = 0.0                  
p_soil = -0.001
p_base_ini = -4.0
p_leaf_ini = -5.0

theta_base_0 = p_to_theta(p_base_ini, pi_0_stem, pi_tlp_stem, theta_r_stem, theta_ft_stem, epsilon_stem, S_stem)
theta_leaf_0 = p_to_theta(p_leaf_ini, pi_0_leaf, pi_tlp_leaf, theta_r_leaf, theta_ft_leaf, epsilon_leaf, S_leaf)
y1_0 = float(theta_base_0*size_reservoir_stem_moles)
y2_0 = float(theta_leaf_0*size_reservoir_leaf_moles)
y0 = [y1_0; y2_0]

############### Simulation length ###
tend = 60*60*2
tspan = (0.0,tend)
dt = 1
alg = Euler()
alg_name = "Euler"

################ Solve the problem ##
prob = ODEProblem(roots,y0,tspan,paramset)
sol = solve(prob,alg,adaptive=false,dt=dt)

y_1 = reduce(hcat,sol.u)[1,:]
y_2 = reduce(hcat,sol.u)[2,:]

plot(sol.t, y_1, label="stem", xaxis="t [s]", yaxis="water content [mol]", dpi=500)
plot!(sol.t, y_2, label="leaves")
savefig("test1_water_content_moles.png") 

# Convert soln to volumetric water content and plot
y_theta_1 = y_1/size_reservoir_stem_moles
y_theta_2 = y_2/size_reservoir_leaf_moles

plot(sol.t, y_theta_1, label="stem", xaxis="t [s]", yaxis="relative water content [mol/mol]",dpi=500)
plot!(sol.t,y_theta_2, label="leaves",dpi=500)
savefig("test1_relative_water_content.png") 

# Compute pressure and flow rates from soln and plot
p_base = theta_to_p(y_theta_1, pi_0_stem, theta_r_stem, theta_tlp_stem, theta_ft_stem, epsilon_stem, S_stem)       
p_leaf = theta_to_p(y_theta_2, pi_0_leaf, theta_r_leaf, theta_tlp_leaf, theta_ft_leaf, epsilon_leaf, S_leaf)

p_dos = ones(length(sol.t),1)
p_ups = ones(length(sol.t),1)
sign_flow = ones(length(sol.t),1)
sign_p = ones(length(sol.t),1)

for i in 1:length(sol.t)
    if p_base[i] < p_soil 
        p_dos[i] = p_base[i]
        p_ups[i] = p_soil
        sign_p[i] = 1
        pdif = sign_p[i]*(p_ups[i] - p_dos[i])
        grav = -1*rhog_MPa* h_root/2
        sumf = pdif + grav
        if sumf >= 0
            sign_flow[i] = 1
        else
            sign_flow[i] = -1
        end
    else
        p_dos[i] = p_soil
        p_ups[i] = p_base[i]
        sign_p[i] = -1
        pdif = sign_p[i]*(p_ups[i] - p_dos[i])
        grav = -1*rhog_MPa* h_root/2
        sumf = pdif + grav
        sign_flow[i] = -1
    end
end

flow_in_base_approx = K_max_root_total_moles*vc_integral_approx(p_dos, p_ups, a_root, b_root).*(p_ups .- p_dos .- rhog_MPa.* sign_p.* h_root./2) ./ (p_ups .- p_dos)
flow_in_base = sign_flow.*vc_integral(p_dos, p_ups, sign_p.*h_root/2, flow_in_base_approx, K_max_root_total_moles, rhog_MPa, a_root, b_root)   

p_dos = ones(length(sol.t),1)
p_ups = ones(length(sol.t),1)
sign_flow = ones(length(sol.t),1)
sign_p = ones(length(sol.t),1)

for i in 1:length(sol.t)
    if p_leaf[i] < p_base[i]
        p_dos[i] = p_leaf[i]
        p_ups[i] = p_base[i]
        sign_p[i] = 1
        pdif = sign_p[i]*(p_ups[i] - p_dos[i])
        grav = -1*rhog_MPa* h_stem/2
        sumf = pdif + grav
        if sumf >= 0
            sign_flow[i] = 1
        else
            sign_flow[i] = -1
        end
    else
        p_dos[i] = p_base[i]
        p_ups[i] = p_leaf[i]
        sign_p[i] = -1
        pdif = sign_p[i]*(p_ups[i] - p_dos[i])
        grav = -1*rhog_MPa* h_stem/2
        sumf = pdif + grav
        sign_flow[i] = -1
    end
end

flow_out_base_approx = K_max_stem_moles*vc_integral_approx(p_dos, p_ups, a_stem, b_stem).*(p_ups .- p_dos .- rhog_MPa.* sign_p.*h_stem/2) ./ (p_ups .- p_dos)
flow_out_base = sign_flow.*vc_integral(p_dos, p_ups, sign_p.*h_stem/2, flow_out_base_approx, K_max_stem_moles, rhog_MPa, a_stem, b_stem)  

T = zeros(length(sol.t),1)
p_soil = p_soil.*ones(length(sol.t),1)

# Plot pressure in stem and leaves as function of time [MPa]
plot(sol.t,p_base,linewidth=2,xaxis="time [s]",yaxis="pressure [MPa]",label="stem", dpi=500)
plot!(sol.t,p_leaf,linewidth=2,label="leaves", dpi=500)
plot!(sol.t,p_soil,linewidth=2,label="soil", dpi=500)
savefig("test1_pressure.png") 

# Plot flow as a function of time [mol s-1]
plot(sol.t,flow_in_base,linewidth=2,xaxis="time [s]",yaxis="flow [mol s-1]",label="flow stem",legend=:topright,dpi=500)
plot!(sol.t,flow_out_base,linewidth=2,label="flow leaves",dpi=500)
plot!(sol.t,T,linewidth=2,label="transpiration boundary condition",dpi=500)
savefig("test1_flow.png") 




