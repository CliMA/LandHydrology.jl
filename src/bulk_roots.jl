import XLSX
using OrdinaryDiffEq
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
y0 = [y1_0; y2_0]
tend = 15
tspan = (0,tend)
dt = 1.0

# dt_vec = [0.005, 0.006, 0.007, 0.008, 0.01, 0.05, 0.1, 1, 2];

########## Solver ###########

function roots(dy,y,paramset,t)
    p_base = theta_to_p(y[1]/size_reservoir_stem_moles, pi_0_stem, theta_r_stem, theta_tlp_stem, theta_ft_stem, epsilon_stem, S_stem)      
    p_leaf = theta_to_p(y[2]/size_reservoir_leaf_moles, pi_0_leaf, theta_r_leaf, theta_tlp_leaf, theta_ft_leaf, epsilon_leaf, S_leaf)
      
    flow_in_base_approx = K_max_root_total_moles*vc_integral_approx(p_base, p_soil, a_root, b_root)*(p_soil - p_base - rhog_MPa* h_root/2) / (p_soil - p_base)
    flow_in_base = vc_integral(p_base, p_soil, h_root/2, flow_in_base_approx,K_max_root_total_moles, rhog_MPa, a_root, b_root)   
    flow_out_base_approx = K_max_stem_moles*vc_integral_approx(p_leaf, p_base, a_stem, b_stem)*(p_base - p_leaf - rhog_MPa* h_stem/2) / (p_base - p_leaf)
    flow_out_base = vc_integral(p_leaf, p_base, h_stem/2, flow_out_base_approx,K_max_stem_moles, rhog_MPa, a_stem, b_stem)  
    
    if t < 5 
        T = T_ini;
    elseif t >= 5 && t < 10
        T = 1*(t-5)+T_ini
    else t >= 10
        T = 1*5+T_ini
    end
    dy[1] = flow_in_base - flow_out_base
    dy[2] = flow_out_base - T
end

prob = ODEProblem(roots,y0,tspan,paramset)
alg = Euler()
sol = solve(prob,alg,adaptive=false,dt=1)

y_1 = reduce(hcat,sol.u)[1,:]
y_2 = reduce(hcat,sol.u)[2,:]


matlab_y_1 = 
[10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6085656517
10895.6073742940
10895.6041493249
10895.5978864369
10895.5875898104
10895.5724426541
10895.5523160682
10895.5272521528
10895.4972926829
10895.4624791110]

y_theta_1 = y_1/size_reservoir_stem_moles
y_theta_2 = y_2/size_reservoir_leaf_moles

plot(sol.t,y_theta_1)
plot!(sol.t,y_theta_2)

plot(sol.t, y_1)
plot!(sol.t, matlab_y_1)

#=
p_base = theta_to_p(y_theta_1, pi_0_stem, theta_r_stem, theta_tlp_stem, theta_ft_stem, epsilon_stem, S_stem)       
p_leaf = theta_to_p(y_theta_2, pi_0_leaf, theta_r_leaf, theta_tlp_leaf, theta_ft_leaf, epsilon_leaf, S_leaf)
flow_in_base_approx = K_max_root_total_moles.*vc_integral_approx(p_base, p_soil, a_root, b_root).*(p_soil .- p_base .- rhog_MPa * h_root/2) ./ (p_soil .- p_base)
flow_in_base = vc_integral(p_base, p_soil, h_root/2, flow_in_base_approx,K_max_root_total_moles, rhog_MPa, a_root, b_root)   
flow_out_base_approx = K_max_stem_moles.*vc_integral_approx(p_leaf, p_base, a_stem, b_stem).*(p_base .- p_leaf .- rhog_MPa*h_stem/2) ./ (p_base .- p_leaf)
flow_out_base = vc_integral(p_leaf, p_base, h_stem/2, flow_out_base_approx, K_max_stem_moles, rhog_MPa, a_stem, b_stem)  

T = collect(0.0:Float64(tend))
trange = 1:length(T)
for t in trange
if t < 5 
    T[t] = T_ini;
elseif t >= 5 && t < 10
    T[t] = 1*(t-5)+T_ini
else t >= 10
    T[t] = 1*5+T_ini
end
end

# Plot pressure in stem and leaves as function of time [MPa]
plot(sol.t,p_base,linewidth=2,xaxis="time [s]",yaxis="pressure [MPa]",label="stem")
plot!(sol.t,p_leaf,linewidth=2,label="leaves")

# Plot flow as a function of time [mol s-1]
plot(sol.t,flow_in_base,linewidth=2,xaxis="time [s]",yaxis="flow [mol s-1]",label="flow into stem")
plot!(sol.t,flow_out_base,linewidth=2,label="flow into leaves")
plot!(sol.t,T,linewidth=2,label="transpiration boundary condition")

# Plot error as function of time [mol]
yt = y_1 .+ y_2
ychange = diff(yt)
su = (flow_in_base.-T).*dt
water_cons = ychange.-su[2:end]
tvec=collect(1:length(ychange))
plot(tvec,water_cons,xaxis="time [s]",yaxis="water conservation error per step [mol]", label="julia")
# yaxis("diff(yt)-[(flow in base-T).*dt](2:end) [mol]")

water_cons_matlab = [-8.21565038222616e-15
-8.21565038222616e-15
-8.21565038222616e-15
-8.21565038222616e-15
-8.21565038222616e-15
0.499999248645697
0.499995753433959
0.499989288209396
0.499979878043071
0.499967547822279
-4.69263787996610e-05
-6.15275986266539e-05
-7.60053514987646e-05
-9.03605941262597e-05
-0.000104594268391267]

plot!(tvec,water_cons_matlab,xaxis="time [s]",yaxis="water conservation error per step [mol]",label="matlab")

yt_matlab =[26747.2358918095
26747.2358918095
26747.2358918095
26747.2358918095
26747.2358918095
26747.2358918095
26746.7358920601
26745.2358958122
26742.7359120391
26739.2359556373
26734.7360473766
26729.7362135998
26724.7364834236
26719.7368859681
26714.7374501071
26709.7382044698]

tvec2=collect(1:length(yt))
plot(tvec2,yt)
plot!(tvec2, yt_matlab, xaxis("t"), yaxis("total water content moles"))
=#

#=
matlab_soln = 
[10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087191966
10895.6086000608
10895.6082775635
10895.6076512734
10895.6066216074
10895.6051068844
10895.6030942129
10895.6005878012
10895.5975918253
10895.5941104289
10895.5901477237
10895.5857077896
10895.5807946750
10895.5754123971
10895.5695649419
10895.5632562648
10895.5564902903
10895.5492709131
10895.5416019972
10895.5334873773
10895.5249308579
10895.5159362147
10895.5065071936
10895.4966475119
10895.4863608582
10895.4756508924
10895.4645212461
10895.4529755229
10895.4410172985
10895.4286501210
10895.4158775109
10895.4027029617
10895.3891299395
10895.3751618840
10895.3608022080
10895.3460542981
10895.3309215144
10895.3154071913
10895.2995146372
10895.2832471350
10895.2666079420
10895.2496002906
10895.2322273878
10895.2144924160
10895.1963985329
10895.1779488718
10895.1591465416
10895.1399946272
10895.1204961896
10895.1006542660
10895.0804718703
10895.0599519927
10895.0390976006
10895.0179116380
10894.9963970264
10894.9745566645
10894.9523934285
10894.9299101725
10894.9071097283
10894.8839949057
10894.8605684928
10894.8368332561
10894.8127919406
10894.7884472698
10894.7638019465
10894.7388586522
10894.7136200476
10894.6880887730
10894.6622674478
10894.6361586715
10894.6097650232
10894.5830890618
10894.5561333268
10894.5289003375
10894.5013925938
10894.4736125764
10894.4455627464
10894.4172455460
10894.3886633983
10894.3598187076
10894.3307138595
10894.3013512211
10894.2717331412
10894.2418619500
10894.2117399600
10894.1813694655
10894.1507527428
10894.1198920508
10894.0887896307
10894.0574477063
10894.0258684839
10893.9940541529
10893.9620068856
10893.9297288372
10893.8972221463]
plot!(0:100, matlab_soln)
l2norm = sqrt(sum((matlab_soln .- sol.u)^2))
plot(sol,linewidth=2,xaxis="t",label=["θ [rad]" "ω [rad/s]"],layout=(2,1))
=#

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
