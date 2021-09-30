import XLSX
using OrdinaryDiffEq
using Plots
using SymPy
using LsqFit
using NLsolve
xf = XLSX.readxlsx("parameters_roots.xlsx")
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

# have a default grass, default plant, default tree ; different methods for each type once we start generalizing
# use similar structure as for soil, think about hierarchy and regrouping

#=
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

paramset = [a_root a_stem b_root b_stem K_max_root_total_moles K_max_stem_moles rho_water rhog_MPa h_root h_stem pi_0_leaf pi_0_stem epsilon_leaf epsilon_stem theta_r_leaf theta_r_stem WD theta_tlp_leaf theta_tlp_stem theta_ft_stem size_reservoir_stem_moles size_reservoir_leaf_moles theta_ft_leaf volume_mole_water S_stem S_leaf porosity_stem porosity_leaf];

########## Functions #########
 
    function vc_integral_approx(p_dos, p_ups, a, b)
    #     if p_ups > 0 && p_dos > 0
    #         error('positive pressure value in approx integral')
    #     end
        t_dos = log(a .* exp(b.*p_dos) + 1);
        t_ups = log(a .* exp(b.*p_ups) + 1);
        flux_approx = (a+1) / (a*b) .* (t_ups - t_dos);
        return flux_approx
    end
    
    function vc_integral(p_dos, p_ups, h, flux_approx, Kmax, rhog_MPa, a, b)
    #     if p_ups <= 0 && p_dos <= 0
    #         error('positive pressure value in approx integral')
    #     end
        krghe = Kmax * rhog_MPa * h * (a+1) / a + flux_approx;
        lower = b * krghe;
        multi = Kmax * (a+1) / a .* flux_approx;
        upper_dos = log(a * krghe .* exp(b.*p_dos) + flux_approx);
        upper_ups = log(a * krghe .* exp(b.*p_ups) + flux_approx);
        flux = multi .* (upper_ups - upper_dos) ./ lower;
        return flux
    end
     
    function theta_to_p(theta, pi_0, theta_r, theta_tlp, theta_ft, epsilon, S_s)      
        p = (theta-1)*5;  
        return p
    end
    
    function p_to_theta(p, pi_0, pi_tlp, theta_r, theta_ft, epsilon, S_s)
        theta = p/5 + 1; 
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

########## Initial values ###########
T_ini = 0.01/mass_mole_water # moles per s, using value at noon fig 9j) Christoffersen                    
p_soil = 0 # MPa want it to be wet and wetter than rest of plant, so close to 0

# Set system to equilibrium state by setting LHS of both odes to 0
function fbase!(F_base,x)
    flow_in_base_approx = K_max_root_total_moles*vc_integral_approx(x[1], p_soil, a_root, b_root)*(p_soil - x[1] - rhog_MPa* h_root/2) / (p_soil - x[1])
    flow_in_base = vc_integral(x[1], p_soil, h_root/2, flow_in_base_approx,K_max_root_total_moles, rhog_MPa, a_root, b_root)
    F_base[1] = flow_in_base - T_ini   
end
solnbase = nlsolve(fbase!, [-1.0])
p_base_ini = solnbase.zero[1]

function fleaf!(F_leaf,y)
    flow_out_base_approx = K_max_stem_moles*vc_integral_approx(y[1], p_base_ini, a_stem, b_stem)*(p_base_ini - y[1] - rhog_MPa* h_stem/2) / (p_base_ini - y[1])
    flow_out_base = vc_integral(y[1], p_base_ini, h_stem/2, flow_out_base_approx, K_max_stem_moles, rhog_MPa, a_stem, b_stem)
    F_leaf[1] = flow_out_base - T_ini 
end
solnleaf= nlsolve(fleaf!, [-1.0])
p_leaf_ini = solnleaf.zero[1]

theta_base_0 = p_to_theta(p_base_ini, pi_0_stem, pi_tlp_stem, theta_r_stem, theta_ft_stem, epsilon_stem, S_stem)
theta_leaf_0 = p_to_theta(p_leaf_ini, pi_0_leaf, pi_tlp_leaf, theta_r_leaf, theta_ft_leaf, epsilon_leaf, S_leaf)
y1_0 = float(theta_base_0*size_reservoir_stem_moles)
y2_0 = float(theta_leaf_0*size_reservoir_leaf_moles)
y0 = [y1_0; y2_0; T_ini; T_ini; T_ini]
tspan = (0.0,100.0)
dt = 1.0

# dt_vec = [0.005, 0.006, 0.007, 0.008, 0.01, 0.05, 0.1, 1, 2];
# dt = 1; 

########## Solver ###########

function roots(dy,y,paramset,t)

    #y_theta_1 = y[1]/size_reservoir_stem_moles
    #y_theta_2 = y[2]/size_reservoir_leaf_moles

    p_base = theta_to_p(y[1]/size_reservoir_stem_moles, pi_0_stem, theta_r_stem, theta_tlp_stem, theta_ft_stem, epsilon_stem, S_stem)   
    
    p_leaf = theta_to_p(y[2]/size_reservoir_leaf_moles, pi_0_leaf, theta_r_leaf, theta_tlp_leaf, theta_ft_leaf, epsilon_leaf, S_leaf)
      
    flow_in_base_approx = K_max_root_total_moles*vc_integral_approx(p_base, p_soil, a_root, b_root)*(p_soil - p_base - rhog_MPa* h_root/2) / (p_soil - p_base)
    flow_in_base = vc_integral(p_base, p_soil, h_root/2, flow_in_base_approx,K_max_root_total_moles, rhog_MPa, a_root, b_root)   
    flow_out_base_approx = K_max_stem_moles*vc_integral_approx(p_leaf, p_base, a_stem, b_stem)*(p_base - p_leaf - rhog_MPa* h_stem/2) / (p_base - p_leaf)
    flow_out_base = vc_integral(p_leaf, p_base, h_stem/2, flow_out_base_approx,K_max_stem_moles, rhog_MPa, a_stem, b_stem)  
    
    if t < 5 
        T = T_ini;
    elseif t >= 5 && t < 10
        T = 0.1*(t-5)+T_ini
    else t >= 10
        T = 0.1*5+T_ini
    end
    dy[1] = flow_in_base - flow_out_base
    dy[2] = flow_out_base - T
end

prob = ODEProblem(roots,y0,tspan,p)
alg = Euler()
sol = solve(prob,alg,adaptive=false,dt=1.0)
plot(sol)
#plot(sol,linewidth=2,xaxis="t",label=["θ [rad]" "ω [rad/s]"],layout=(2,1))