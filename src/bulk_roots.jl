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
p_soil = 0.0 # MPa want it to be wet and wetter than rest of plant, so close to 0

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

############### Simulation length and dt ############
tend = 15.0
tspan = (0.0,tend)

# MATLAB soln dt = 1
matlab_y_1 = [10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6077104577
10895.6046417424
10895.5985216132
10895.5883501251
10895.5731358228
10895.5529214613
10895.5277494643
10895.4976619264]

prob = ODEProblem(roots,y0,tspan,paramset)
alg = Euler()
sol = solve(prob,alg,adaptive=false,dt=1)

y_1 = reduce(hcat,sol.u)[1,:]
y_2 = reduce(hcat,sol.u)[2,:]

plot(sol.t, y_1, label="julia, dt = 1s", xaxis="t [s]", yaxis="water content stem [mol]", title="MATLAB vs Julia")
plot!(sol.t, matlab_y_1, seriestype = :scatter, label="matlab, dt = 1s")

# MATLAB soln dt = 0.5
matlab_y_1 = [10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6086080322
10895.6082238999
10895.6074567186
10895.6061798850
10895.6042673326
10895.6015935293
10895.5980334756
10895.5934627023
10895.5877572691
10895.5807937620
10895.5725775153
10895.5631138425
10895.5524080364
10895.5404653696
10895.5272910939
10895.5128904409
10895.4972686219
10895.4804308281]

prob = ODEProblem(roots,y0,tspan,paramset)
alg = Euler()
sol = solve(prob,alg,adaptive=false,dt=0.5)

y_1 = reduce(hcat,sol.u)[1,:]
y_2 = reduce(hcat,sol.u)[2,:]

plot!(sol.t, y_1, label="julia, dt = 0.5s")
plot!(sol.t, matlab_y_1, seriestype = :scatter, label="matlab, dt = 0.5s")

# MATLAB soln dt = 0.1
matlab_y_1 = [10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087362571
10895.6087352313
10895.6087321548
10895.6087260035
10895.6087157542
10895.6087003845
10895.6086788731
10895.6086501992
10895.6086133433
10895.6085672863
10895.6085110104
10895.6084434982
10895.6083637337
10895.6082707012
10895.6081633864
10895.6080407754
10895.6079018554
10895.6077456144
10895.6075710413
10895.6073771258
10895.6071628586
10895.6069272310
10895.6066692353
10895.6063878648
10895.6060821133
10895.6057509757
10895.6053934479
10895.6050085262
10895.6045952082
10895.6041524921
10895.6036793770
10895.6031748628
10895.6026379505
10895.6020676416
10895.6014629386
10895.6008228449
10895.6001463648
10895.5994325032
10895.5986802660
10895.5978886600
10895.5970566927
10895.5961833727
10895.5952677091
10895.5943087120
10895.5933053924
10895.5922567621
10895.5911618337
10895.5900196207
10895.5888291374
10895.5875893989
10895.5862994213
10895.5849592470
10895.5835689187
10895.5821284790
10895.5806379703
10895.5790974352
10895.5775069160
10895.5758664552
10895.5741760951
10895.5724358782
10895.5706458467
10895.5688060428
10895.5669165089
10895.5649772872
10895.5629884198
10895.5609499488
10895.5588619165
10895.5567243649
10895.5545373359
10895.5523008718
10895.5500150143
10895.5476798056
10895.5452952875
10895.5428615018
10895.5403784905
10895.5378462954
10895.5352649583
10895.5326345209
10895.5299550249
10895.5272265120
10895.5244490240
10895.5216226024
10895.5187472888
10895.5158231248
10895.5128501520
10895.5098284117
10895.5067579456
10895.5036387949
10895.5004710012
10895.4972546058
10895.4939896500
10895.4906761752
10895.4873142226
10895.4839038335
10895.4804450490
10895.4769379104
10895.4733824588
10895.4697787352
10895.4661267809]

dt = 0.1
prob = ODEProblem(roots,y0,tspan,paramset)
alg = Euler()
sol = solve(prob,alg,adaptive=false,dt=dt)

y_1 = reduce(hcat,sol.u)[1,:]
y_2 = reduce(hcat,sol.u)[2,:]

plot!(sol.t, y_1, label="julia, dt = 0.1s")
plot!(sol.t, matlab_y_1, seriestype = :scatter, label="matlab, dt = 0.1s")

# Convert soln to volumetric water content and plot
y_theta_1 = y_1/size_reservoir_stem_moles
y_theta_2 = y_2/size_reservoir_leaf_moles

plot(sol.t,y_theta_1)
plot!(sol.t,y_theta_2)

# Compute pressure and flow rates from soln and plot
p_base = theta_to_p(y_theta_1, pi_0_stem, theta_r_stem, theta_tlp_stem, theta_ft_stem, epsilon_stem, S_stem)       
p_leaf = theta_to_p(y_theta_2, pi_0_leaf, theta_r_leaf, theta_tlp_leaf, theta_ft_leaf, epsilon_leaf, S_leaf)
flow_in_base_approx = K_max_root_total_moles.*vc_integral_approx(p_base, p_soil, a_root, b_root).*(p_soil .- p_base .- rhog_MPa * h_root/2) ./ (p_soil .- p_base)
flow_in_base = vc_integral(p_base, p_soil, h_root/2, flow_in_base_approx,K_max_root_total_moles, rhog_MPa, a_root, b_root)   
flow_out_base_approx = K_max_stem_moles.*vc_integral_approx(p_leaf, p_base, a_stem, b_stem).*(p_base .- p_leaf .- rhog_MPa*h_stem/2) ./ (p_base .- p_leaf)
flow_out_base = vc_integral(p_leaf, p_base, h_stem/2, flow_out_base_approx, K_max_stem_moles, rhog_MPa, a_stem, b_stem)  

T = collect(0.0:dt:Float64(tend))
tlength = 1:length(T)
for i in tlength
if T[i] < 5
    T[i] = T_ini;
elseif T[i] >= 5 && T[i] < 10
    T[i] = 1*(T[i]-5)+T_ini
else T[i] >= 10
    T[i] = 1*5+T_ini
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
water_cons_matlab = [-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
-8.21565038222616e-16
0.00999999999999923
0.0100000000015999
0.00999999939711934
0.00999999819580194
0.00999999638557426
0.00999999398397851
0.00999999097268399
0.00999998736297104
0.00999998315935412
0.00999997835594065
0.00999997295098098
0.00999996695051234
0.00999996035016013
0.00999995315333588
0.00999994536395293
0.00999993697551782
0.00999992799295690
0.00999991841442341
0.00999990823857533
0.00999989747548807
0.00999988611755001
0.00999987416129269
0.00999986161829877
0.00999984847883009
0.00999983475183430
0.00999982043585307
0.00999980552264729
0.00999979002631152
0.00999977393960710
0.00999975726307678
0.00999973999412340
0.00999972214156131
0.00999970370379022
0.00999968467607237
0.00999966506144073
0.00999964486343169
0.00999962407515609
0.00999960270842015
0.00999958075733354
0.00999955821286519
0.00999953509558060
0.00999951139380118
0.00999948710726034
0.00999946223618164
0.00999943678856075
0.00999941075833510
0.00999938414721147
0.00999935695739196
0.00999932918429408
0.00999930083837985
-7.28091224366789e-07
-7.57591035571714e-07
-7.87074179053882e-07
-8.16530206004806e-07
-8.45959593043233e-07
-8.75366480412865e-07
-9.04744113183753e-07
-9.34096665794471e-07
-9.63431974199125e-07
-9.92732425153520e-07
-1.02201680574066e-06
-1.05127118099624e-06
-1.08050346181443e-06
-1.10971066358267e-06
-1.13888982572474e-06
-1.16804527716674e-06
-1.19717772961714e-06
-1.22628427834393e-06
-1.25536931094850e-06
-1.28442632207904e-06
-1.31345973741759e-06
-1.34246672850447e-06
-1.37144812356604e-06
-1.40040840379552e-06
-1.42934479974155e-06
-1.45825455799553e-06
-1.48713858072469e-06
-1.51599779102396e-06
-1.54483676906292e-06
-1.57364556185868e-06
-1.60243606284194e-06
-1.63119471729978e-06
-1.65993709860812e-06
-1.68864968697058e-06
-1.71733717918565e-06
-1.74600792346435e-06
-1.77464118839010e-06
-1.80326354964899e-06
-1.83185431290589e-06
-1.86042190075675e-06
-1.88896384906156e-06
-1.91748497813116e-06
-1.94597922242856e-06
-1.97444780619671e-06
-2.00289561175238e-06
-2.03131662662726e-06
-2.05971577016406e-06
-2.08809070395555e-06
-2.11643911229853e-06
-2.14476596904767e-06]

T_matlab = [0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.555084350617920
0.655084350617920
0.755084350617920
0.855084350617921
0.955084350617920
1.05508435061792
1.15508435061792
1.25508435061792
1.35508435061792
1.45508435061792
1.55508435061792
1.65508435061792
1.75508435061792
1.85508435061792
1.95508435061792
2.05508435061792
2.15508435061792
2.25508435061792
2.35508435061792
2.45508435061792
2.55508435061792
2.65508435061792
2.75508435061792
2.85508435061792
2.95508435061792
3.05508435061792
3.15508435061792
3.25508435061792
3.35508435061792
3.45508435061792
3.55508435061792
3.65508435061792
3.75508435061792
3.85508435061792
3.95508435061792
4.05508435061792
4.15508435061792
4.25508435061792
4.35508435061792
4.45508435061792
4.55508435061792
4.65508435061792
4.75508435061792
4.85508435061792
4.95508435061792
5.05508435061792
5.15508435061792
5.25508435061792
5.35508435061792
5.45508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792
5.55508435061792]

plot(sol.t,T,xaxis="time [s]",yaxis="T [mol s-1]", label="julia")
plot!(sol.t,T_matlab,xaxis="time [s]",yaxis="T [mol s-1]",label="matlab")

yt = y_1 .+ y_2
ychange = diff(yt)
su = (flow_in_base.-T).*dt
water_cons = ychange.-su[2:end]
tvec=collect(1:length(ychange))

plot(tvec,water_cons,xaxis="time [s]",yaxis="water conservation error per step [mol]", label="julia", title="diff(yt)-[(flow in base-T).*dt](2:end) [mol]")
plot!(tvec,water_cons_matlab,seriestype = :scatter, xaxis="time [s]",yaxis="water conservation error per step [mol]",label="matlab")

#=
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
