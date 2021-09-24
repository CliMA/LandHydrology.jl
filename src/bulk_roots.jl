clc;
clear all;
close all;
tic
%% Constants and parameters
T = table2cell(readtable('parameters_roots.xlsx'));
params = readmatrix('parameters_roots.xlsx','Range',['C2:C', num2str(size(T,1)+1)]);
MPa_to_Pa = params(find(strcmp(T(:,2),'MPa_to_Pa')));
rho_water = params(find(strcmp(T(:,2),'rho_water')));
rhog_MPa = params(find(strcmp(T(:,2),'rhog_MPa')));
mass_mole_water = params(find(strcmp(T(:,2),'mass_mole_water'))); 
volume_mole_water = params(find(strcmp(T(:,2),'volume_mole_water'))); 
h_root = params(find(strcmp(T(:,2),'h_root'))); 
h_stem = params(find(strcmp(T(:,2),'h_stem')));
LMA = params(find(strcmp(T(:,2),'LMA')));
WD = params(find(strcmp(T(:,2),'WD')));
K_max_stem = params(find(strcmp(T(:,2),'K_max_stem')));
K_max_root = params(find(strcmp(T(:,2),'K_max_root')));
K_max_root_total = params(find(strcmp(T(:,2),'K_max_root_total')));
K_max_stem_moles = params(find(strcmp(T(:,2),'K_max_stem_moles')));
K_max_root_moles = params(find(strcmp(T(:,2),'K_max_root_moles')));
K_max_root_total_moles = params(find(strcmp(T(:,2),'K_max_root_total_moles')));
porosity_leaf = params(find(strcmp(T(:,2),'porosity_leaf')));
porosity_stem = params(find(strcmp(T(:,2),'porosity_stem')));
size_reservoir_stem = params(find(strcmp(T(:,2),'size_reservoir_stem')));
size_reservoir_leaf = params(find(strcmp(T(:,2),'size_reservoir_leaf')));
size_reservoir_stem_moles = params(find(strcmp(T(:,2),'size_reservoir_stem_moles')));
size_reservoir_leaf_moles = params(find(strcmp(T(:,2),'size_reservoir_leaf_moles')));
C_stem = params(find(strcmp(T(:,2),'C_stem')));
C_leaf = params(find(strcmp(T(:,2),'C_leaf')));
S_stem = params(find(strcmp(T(:,2),'S_stem')));
S_leaf = params(find(strcmp(T(:,2),'S_leaf')));

epsilon_leaf = params(find(strcmp(T(:,2),'epsilon_leaf')));
epsilon_stem = params(find(strcmp(T(:,2),'epsilon_stem')));
pi_0_leaf = params(find(strcmp(T(:,2),'pi_o_leaf')));
pi_0_stem = params(find(strcmp(T(:,2),'pi_o_stem')));
pi_tlp_leaf = params(find(strcmp(T(:,2),'pi_tlp_leaf')));
pi_tlp_stem = params(find(strcmp(T(:,2),'pi_tlp_stem')));
theta_r_leaf = params(find(strcmp(T(:,2),'theta_r_leaf')));
theta_r_stem = params(find(strcmp(T(:,2),'theta_r_stem')));
theta_tlp_leaf = params(find(strcmp(T(:,2),'theta_tlp_leaf')));
theta_tlp_stem = params(find(strcmp(T(:,2),'theta_tlp_stem')));
theta_ft_stem = params(find(strcmp(T(:,2),'theta_ft_stem')));
theta_ft_leaf = params(find(strcmp(T(:,2),'theta_ft_leaf')));
f_cap = params(find(strcmp(T(:,2),'f_cap')));
P_50_stem = params(find(strcmp(T(:,2),'P_50_stem')));
P_50_root = params(find(strcmp(T(:,2),'P_50_root')));
a_x_stem = params(find(strcmp(T(:,2),'a_x_stem')));
a_x_root = params(find(strcmp(T(:,2),'a_x_root')));

%% Find which theta_to_p conversion is physical for 2nd regime
% syms pi_0 epsilon theta_r theta_ft p theta
% p_sol = -abs(pi_0)*(theta_ft-theta_r)/(theta-theta_r);
% p_p = abs(pi_0)-epsilon*(theta_ft-theta)/(theta_ft-theta_r);
% eqn = p_sol + p_p == p;
% theta_as_fnc_of_p = solve(eqn, theta); 
 
% pi_0 = pi_0_leaf;
% theta_ft = theta_ft_leaf;
% theta_tlp = theta_tlp_leaf;
% theta_r = theta_r_leaf;
% epsilon = epsilon_leaf;
% p_tlp = p_tlp_leaf;
% p_ft = p_ft_leaf;
% p = linspace(p_tlp,p_ft,50);
 
% figure(1)
% theta_1 = (epsilon*theta_r + epsilon*theta_ft - p.*theta_r + p.*theta_ft + theta_r*abs(pi_0) - theta_ft*abs(pi_0) + theta_r*(abs(pi_0).^2 + 2*epsilon.*p + 2*epsilon*abs(pi_0) - 2.*p.*abs(pi_0) + epsilon^2 + p.^2).^(1/2) - theta_ft*(abs(pi_0)^2 + 2*epsilon.*p + 2*epsilon*abs(pi_0) - 2.*p*abs(pi_0) + epsilon^2 + p.^2).^(1/2))/(2*epsilon);
% theta_2 = (epsilon*theta_r + epsilon*theta_ft - p.*theta_r + p.*theta_ft + theta_r*abs(pi_0) - theta_ft*abs(pi_0) - theta_r*(abs(pi_0).^2 + 2*epsilon.*p + 2*epsilon*abs(pi_0) - 2.*p.*abs(pi_0) + epsilon^2 + p.^2).^(1/2) + theta_ft*(abs(pi_0)^2 + 2*epsilon.*p + 2*epsilon*abs(pi_0) - 2.*p*abs(pi_0) + epsilon^2 + p.^2).^(1/2))/(2*epsilon);
% plot(p, theta_1)
% hold on
% plot(p, theta_2)
% xline(p_tlp, '-g')
% xline(p_ft, '-b')
% legend('theta_1', 'theta_2')
% xlabel('pressure [MPa]')
% ylabel('vwc []')
% legend('soln 1','soln 2','p_{tlp}','p_{ft}')
% axis([-3 1 0.2 1])
% % result : theta 2 is unphysical

%% Plot PV curves for stem. P to theta, and theta to P
% Theta to P
% figure(2)
% theta = (theta_r_stem+ 0.01):0.001:2.5;
% p = nan(1,length(theta));
% for i = 1:length(theta)
%     p(i) = theta_to_p(theta(i), pi_0_stem, theta_r_stem, theta_tlp_stem, theta_ft_stem,...
%         epsilon_stem, S_stem);
% end
% plot(theta, p, '.k')
% hold on
% xline(theta_r_stem, '-r');
% xline(theta_tlp_stem, '-g')
% xline(theta_ft_stem, '-b')
% xlabel('relative water content []')
% ylabel('pressure [MPa]')
% legend('p',['theta_{r_{stem}}=' num2str(theta_r_stem)],...
%     ['theta_{tlp_{stem}}=' num2str(theta_tlp_stem)],['theta_{ft_{stem}}='...
%     num2str(theta_ft_stem)],'location','best')
% axis([0 1.1 -10 2])
% 
% % P to theta
% figure(3)
% theta_from_p = nan(1, length(p));
% for i = 1:length(p)
%     theta_from_p(i) = p_to_theta(p(i), pi_0_stem, ...
%         pi_tlp_stem, theta_r_stem, theta_ft_stem, epsilon_stem, S_stem);
% end
% plot(p, theta_from_p, '.k')
% hold on
% xline(pi_tlp_stem, '-g')
% xline(0, '-b')
% ylabel('relative water content')
% xlabel('pressure [MPa]')
% legend('theta',['pi_{tlp_{stem}}=' num2str(pi_tlp_stem)],['pi_{ft_{stem}}=' num2str(0)],'location','best')
% axis([-10 1 0 1.1])
%  
% % Check the conversion, i.e. that theta = theta_from_p
% % check=all(theta==theta_from_p)==1; 
% figure(4)
% plot(p, theta_from_p, 'or')
% hold on
% plot(p, theta, '-k')
% xlabel('pressure [MPa]')
% ylabel('vwc []')
% title('Stem PV curve')
% legend('theta from p','original theta','location','best')
% axis([-10 1 0 1.1])
%   
%% Plot PV curves for leaf. P to theta, and theta to P
% % % Leaf
% % Theta to P
% figure(5)
% theta = (theta_r_leaf+ 0.01):0.001:2.5;
% p = nan(1,length(theta));
% for i = 1:length(theta)
%     p(i) = theta_to_p(theta(i), pi_0_leaf, theta_r_leaf, theta_tlp_leaf, theta_ft_leaf,...
%         epsilon_leaf, S_leaf);
% end
% plot(theta, p, '.k')
% hold on
% xline(theta_r_leaf, '-r');
% xline(theta_tlp_leaf, '-g')
% xline(theta_ft_leaf, '-b')
% xlabel('relative water content []')
% ylabel('pressure [MPa]')
% legend('p',['theta_{r_{leaf}}=' num2str(theta_r_leaf)],...
%     ['theta_{tlp_{leaf}}=' num2str(theta_tlp_leaf)],['theta_{ft_{leaf}}='...
%     num2str(theta_ft_leaf)],'location','best')
% axis([0 1.1 -10 2])
% 
% % P to theta
% figure(6)
% theta_from_p = nan(1, length(p));
% for i = 1:length(p)
%     theta_from_p(i) = p_to_theta(p(i), pi_0_leaf, ...
%         pi_tlp_leaf, theta_r_leaf, theta_ft_leaf, epsilon_leaf, S_leaf);
% end
% plot(p, theta_from_p, '.k')
% hold on
% xline(pi_tlp_leaf, '-g')
% xline(0, '-b')
% ylabel('relative water content')
% xlabel('pressure [MPa]')
% legend('theta',['pi_{tlp_{leaf}}=' num2str(pi_tlp_leaf)],['pi_{ft_{leaf}}=' num2str(0)],'location','best')
% axis([-10 1 0 1.1])
%  
% % Check the conversion, i.e. that theta = theta_from_p
% % check=all(theta==theta_from_p)==1; 
% figure(7)
% plot(p, theta_from_p, 'or')
% hold on
% plot(p, theta, '-k')
% xlabel('pressure [MPa]')
% ylabel('vwc []')
% title('Leaf PV curve')
% legend('theta from p','original theta','location','best')
% axis([-10 1 0 1.1])

%% Plot k(P) for stem

% Testing Christoffersen's functions for psi_50 and a_x
% figure(8)
% WD=linspace(0,1,100);
% p = linspace(-10,0,1000);
% plot(WD,(-((3.57.*WD).^1.73))-1.09) % Christoffersen's psi_50(WD)
% axis([0 1 -10 1])
% plot(p,(54.4.*(-p).^(-1.17))/100) % Christoffersen's a_x(WD)
% axis([-5 1 0 100])

% figure(9)
P = linspace(-15,0,1000);
k_stem_power = (1 + (P./P_50_stem).^(a_x_stem)).^-1;
logisticEqn = '(a_stem+1)/a_stem*(1-1/(1+a_stem*exp(b_stem*x)))';
startPoints = [1 1];
fstem = fit(P',k_stem_power', logisticEqn,'Start', startPoints); 
stem_coeffs = coeffvalues(fstem);
a_stem = stem_coeffs(1);
b_stem = stem_coeffs(2);
% plot(P,k_stem_power,'-b')
% hold on
% plot(P,(a_stem+1)/a_stem*(1-1./(1+a_stem.*exp(b_stem.*P))),'-r')

%% Plot k(P) for root
k_root_power = (1 + (P./P_50_root).^(a_x_root)).^-1;
logisticEqn = '(a_root+1)/a_root*(1-1/(1+a_root*exp(b_root*x)))';
startPoints = [1 1];
froot = fit(P',k_root_power', logisticEqn,'Start', startPoints); 
root_coeffs = coeffvalues(froot);
a_root = root_coeffs(1);
b_root = root_coeffs(2);
% plot(P,k_root_power,'-k')
% plot(P,(a_root+1)/a_root*(1-1./(1+a_root.*exp(b_root.*P))),'-g')
% % Check that function stays at 1 above P=0
% % P = linspace(-15,5,1000);
% % plot(P,(a_root+1)/a_root*(1-1./(1+a_root.*exp(b_root.*P))),'-r')
% axis([-15 0 0 1])
% xlabel('P [MPa]')
% ylabel('k(P)')
% legend('power vc for stem','fitted logistic vc for stem','power vc for root','fitted logistic vc for root')
% title('Stem and root vulnerability curves')

%% Make sure boundaries of integral are correct. Pdos = pressure downstream, Pups = pressure upstream (aka, tip of roots for us)
% p_soil_vec = [0, -0.5, -1.0, -1.5, -2.5];
% p_base_init = linspace(0, -5, 100);
% figure(20)
% for j = 1:length(p_soil_vec)
%     p_soil=p_soil_vec(j)
%     j
%     for i=p_base_init(1:end)  
%         flow_in_base_approx = K_max_root_total_moles*vc_integral_approx(i, p_soil, a_root, b_root)*...
%             (p_soil - i - rhog_MPa* h_root/2) / (p_soil - i);
%         flow_in_base = vc_integral(i, p_soil, h_root/2, flow_in_base_approx,...
%             K_max_root_total_moles, rhog_MPa, a_root, b_root);
%         plot(i,flow_in_base, '.k')
%         hold on
%         xlabel('p base')
%         ylabel('flow')
%     end
% end

%% Parameters
paramset = [a_root, a_stem, b_root, b_stem, K_max_root_total_moles, K_max_stem_moles,...
    rho_water, rhog_MPa, h_root, h_stem, pi_0_leaf, pi_0_stem, epsilon_leaf,...
    epsilon_stem, theta_r_leaf, theta_r_stem,...
    WD, theta_tlp_leaf, theta_tlp_stem, theta_ft_stem,...
    size_reservoir_stem_moles, size_reservoir_leaf_moles,...
    theta_ft_leaf, volume_mole_water, S_stem, S_leaf, porosity_stem, porosity_leaf];
paramset_names = {'a_root', 'a_stem', 'b_root', 'b_stem', 'K_max_root_total_moles',...
     'K_max_stem_moles','rho_water', 'rhog_MPa', 'h_root', 'h_stem', 'pi_0_leaf', 'pi_0_stem',...
     'epsilon_leaf','epsilon_stem',...
     'theta_r_leaf', 'theta_r_stem', 'WD', 'theta_tlp_leaf',...
     'theta_tlp_stem', 'theta_ft_stem', 'size_reservoir_stem_moles',...
     'size_reservoir_leaf_moles','theta_ft_leaf',...
     'volume_mole_water', 'S_stem', 'S_leaf', 'porosity_stem', 'porosity_leaf'};
 
%% Initial states and length of simulation
T_ini = 0.01/mass_mole_water; % moles per s, using value at noon fig 9j) Christoffersen                    
p_soil = 0; % MPa want it to be wet and wetter than rest of plant, so close to 0
%dt_vec = [0.005, 0.006, 0.007, 0.008, 0.01, 0.05, 0.1, 1, 2];
dt = 0.01; 

% for i=1:length(dt_vec)
%     dt=dt_vec(i)
tspan = 0:dt:15; % s length of simulation 
%% Set system to equilibrium state by setting LHS of both odes to 0
syms p_base_ini
flow_in_base_approx = K_max_root_total_moles*vc_integral_approx(p_base_ini, p_soil, a_root, b_root)*...
(p_soil - p_base_ini - rhog_MPa* h_root/2) / (p_soil - p_base_ini);
flow_in_base = vc_integral(p_base_ini, p_soil, h_root/2, flow_in_base_approx,...
K_max_root_total_moles, rhog_MPa, a_root, b_root);
eqn = flow_in_base == T_ini;
p_base_ini = double(vpasolve(eqn, p_base_ini));
% check that flow in the base = T_ini as supposed to be in equilibrium
flow_in_base = eval(flow_in_base)*mass_mole_water; % great, matches T_ini/mass_mole_water
theta_base_0 = p_to_theta(p_base_ini, pi_0_stem, pi_tlp_stem,...
     theta_r_stem, theta_ft_stem, epsilon_stem, S_stem);
y1_0 = double(theta_base_0*size_reservoir_stem_moles); % moles
 
syms p_leaf_ini
flow_out_base_approx = K_max_stem_moles*vc_integral_approx(p_leaf_ini, p_base_ini, a_stem, b_stem)*...
    (p_base_ini - p_leaf_ini - rhog_MPa* h_stem/2) / (p_base_ini - p_leaf_ini);
flow_out_base = vc_integral(p_leaf_ini, p_base_ini, h_stem/2, flow_out_base_approx,...
    K_max_stem_moles, rhog_MPa, a_stem, b_stem);
eqn2 = flow_out_base == T_ini;
p_leaf_ini = double(vpasolve(eqn2, p_leaf_ini));
% check that flow out of the base = T_ini as supposed to be in equilibrium
flow_out_base = eval(flow_out_base)*mass_mole_water; % great, matches T_ini/mass_mole_water
theta_leaf_0 = p_to_theta(p_leaf_ini, pi_0_leaf, pi_tlp_leaf,...
     theta_r_leaf, theta_ft_leaf, epsilon_leaf, S_leaf);
y2_0 = double(theta_leaf_0*size_reservoir_leaf_moles);
y0 = [y1_0; y2_0];

%% Run model
% % % Adaptive time stepping
% [t,y] = ode45(@(t,y) odefcn(t,y,paramset,paramset_names,T_ini,p_soil), tspan, y0);
%  
% Fixed time stepping
Y = ode5(@(t,y) odefcn(t,y,paramset,paramset_names,T_ini,p_soil), tspan, y0);

% 2 pools
[dydt,p_base,p_leaf,flow_in_base,flow_out_base,T] = cellfun(@(t,y) odefcn(t,y,paramset,paramset_names,T_ini,p_soil),...
   (num2cell(tspan))', num2cell(Y,2), 'uni', 0);
% % 1 pool
% [dydt,p_base,flow_in_base,T] = cellfun(@(t,y) odefcn(t,y,paramset,paramset_names,T_ini,p_soil),...
%      num2cell(t), num2cell(y,1), 'uni', 0);	 
 
%% Plot model results
% p_base = cell2mat(p_base);
% p_leaf = cell2mat(p_leaf);
T = cell2mat(T);
flow_in_base = cell2mat(flow_in_base);
flow_out_base = cell2mat(flow_out_base);

% % Pressure as a function of time
% fig10=figure(10)
% plot(tspan, p_base, '-b')
% hold on
% % yline(p_base_ini) 
% plot(tspan, p_leaf, '-r')
% %plot(t, p_soil*(ones(length(t))), '-k')
% plot(tspan, p_soil*(ones(length(tspan))), '-k')
% legend('p_{stem}','p_{leaf}','p_{soil}','location','best')
% xlabel('time [s]')
% ylabel('pressure [MPa]')
% % axis([0 1.5*t(end) -0.4 0.1])
% % title(['Transpiration bc=', num2str(T_ini), '[mol s-1] ', 'Root pressure bc=', num2str(p_soil), 'MPa ', 'dt=', num2str(dt), 's ', 'sim length=', num2str(tspan(end)), 's'])
% fname=['Results\pressure_as_fnc_of_t_',datestr(now, 'yyyy-MM-dd-HH-mm-ss'),'.jpg'];
% saveas(fig10,fname);
% 
% %% Water in reservoirs as a function of time
% fig11=figure(11)
% plot(tspan,Y(:,1)/size_reservoir_stem_moles, '-b')
% %plot(tspan,Y-y1_0,'-b')
% hold on
% plot(tspan,Y(:,2)/size_reservoir_leaf_moles, '-r')
% legend('theta_{stem}','theta_{leaf}','location','best')
% xlabel('time [s]')
% % ylabel('water content - inital water content [mol]')
% ylabel('water content [mol/mol at saturation]')
% % axis([0 1.5*t(end) -15 0])
% % title('')
% fname=['Results\rwc_as_fnc_of_t_',datestr(now, 'yyyy-MM-dd-HH-mm-ss'),'.jpg'];
% saveas(fig11,fname);

%% Flow as a function of time
% fig12=figure(12)
% plot(tspan, flow_in_base, '.b')
% hold on
% plot(tspan, flow_out_base, '.r')
% plot(tspan, T,'.g')
% xlabel('time [s]')
% ylabel('flow [mol s-1]')
% legend('flow_{in_{stem}}','flow_{out_{stem}}','Transpiration','location','best')
% fname=['Results\flow_as_fnc_of_t_',datestr(now, 'yyyy-MM-dd-HH-mm-ss'),'.jpg'];
% saveas(fig12,fname);

%% Water conservation 
% fig13=figure(13)
% % plot(t, y(:,1)/1e4)
% % hold on
% % plot(t, y(:,2)/1e4)
% y0 = Y(1); %+y(1,2);
% yt = Y(:); %+y(:,2);
% % plot(t,log(abs((yt-y0-cumsum((flow_in_base-T).*dt))./y0))) % log
% % plot(t,(yt-y0-cumsum((flow_in_base-T).*dt))./y0) % relative to initial value
% ychange = diff(yt);
% cusu = (flow_in_base-T).*dt;
% %plot(tspan,ychange - cusu) % absolute
% plot(tspan(2:end),ychange)
% hold on
% plot(tspan(1:end-1),cusu(1:end-1))
% xlabel('time [s]')
% ylabel('yt-y0-cumsum((flow_in_base-T).*dt) [mol]')
% fname=['Results\water_conservation',datestr(now, 'yyyy-MM-dd-HH-mm-ss'),'.jpg'];
% saveas(fig13,fname);
% legend('change in water in stem [mol]','(flow_in_base-T).*dt [mol]','total tree water content (leaf + stem) [mol]')

% fig14=figure(14)
y0 = Y(1,1)+Y(1,2);
yt = Y(:,1)+Y(:,2);
yend = Y(end,1)+Y(end,2);
ychange = yend-y0;
cusu = sum(flow_in_base-T).*dt;
water_cons_end = ychange-cusu;
rel_error = water_cons_end/ychange;

% ychange = (diff(yt));
% su = (flow_in_base-T).*dt;
% water_cons=ychange-su(2:end);
% plot(linspace(1,tspan(end),length(ychange)),water_cons,'.b')
% hold on

% % plot(tspan(2:end),ychange)
% % hold on
% % plot(tspan(1:end-1),cusu(1:end-1))
xlabel('time [s]')
ylabel('diff(yt)-[(flow in base-T).*dt](2:end) [mol]')
fname=['Results\water_conservation_diff_yt',datestr(now, 'yyyy-MM-dd-HH-mm-ss'),'.jpg'];
% saveas(fig14,fname);
% % legend(['dt=',num2str(dt),'s'])
 
fig15=figure(15)
water_cons2=(yt-y0-(cumsum(flow_in_base-T).*dt));
plot(tspan',water_cons2)
% hold on
% end
xlabel('time [s]')
ylabel('water conseration [mol], (yt-y0-(cumsum((flow in base-T).*dt)))')
fname15=['Results\water_conservation_cum',datestr(now, 'yyyy-MM-dd-HH-mm-ss'),'.jpg'];
% legend('0.005s', '0.006', '0.007', '0.008', '0.01', '0.05', '0.1', '1', '2'); 
legend('0.01s')
saveas(fig15,fname15)
% legend(['dt=' num2str(dt)]) %,'0.05','0.1','1','2','location','SE')

% fig16=figure(16)
% N=tspan(end)./dt;
% exp_err=N.*(dt.^5);
% vec_err = exp_err.*ones(1,length(tspan))';
% semilogy(tspan, vec_err)
% % legend('err dt=0.01','err 0.05', 'err 0.1', 'err 1', 'err 2', 'location','SE')à
% % fname=['water_conservation_expected_error',datestr(now, 'yyyy-MM-dd-HH-mm-ss'),'.jpg'];
% % saveas(fig16,fname)

toc

%% Functions
function [dydt,p_base,p_leaf,flow_in_base,flow_out_base,T] = odefcn(t,y,paramset,paramset_names,T_ini,p_soil)
%function [dydt,p_base,flow_in_base,T] = odefcn(t,y,paramset,paramset_names,T_ini,p_soil)
    S_stem = paramset(find(strcmp('S_stem', paramset_names)));
    S_leaf = paramset(find(strcmp('S_leaf', paramset_names)));
    size_reservoir_stem_moles = paramset(find(strcmp('size_reservoir_stem_moles', paramset_names)));
    size_reservoir_leaf_moles = paramset(find(strcmp('size_reservoir_leaf_moles', paramset_names)));
    K_max_stem_moles = paramset(find(strcmp('K_max_stem_moles', paramset_names))); 
    K_max_root_total_moles = paramset(find(strcmp('K_max_root_total_moles', paramset_names)));
    h_stem = paramset(find(strcmp('h_stem', paramset_names)));
    h_root = paramset(find(strcmp('h_root', paramset_names)));
    rhog_MPa = paramset(find(strcmp('rhog_MPa', paramset_names)));
    volume_mol_water = paramset(find(strcmp('volume_mol_water', paramset_names)));
    pi_0_stem = paramset(find(strcmp('pi_0_stem', paramset_names)));
    pi_0_leaf = paramset(find(strcmp('pi_0_leaf', paramset_names)));
    theta_r_stem = paramset(find(strcmp('theta_r_stem', paramset_names)));
    theta_r_leaf = paramset(find(strcmp('theta_r_leaf', paramset_names)));
    theta_tlp_stem = paramset(find(strcmp('theta_tlp_stem', paramset_names)));
    theta_tlp_leaf = paramset(find(strcmp('theta_tlp_leaf', paramset_names)));
    theta_ft_stem = paramset(find(strcmp('theta_ft_stem', paramset_names)));
    theta_ft_leaf = paramset(find(strcmp('theta_ft_leaf', paramset_names)));
    epsilon_stem = paramset(find(strcmp('epsilon_stem', paramset_names)));
    epsilon_leaf = paramset(find(strcmp('epsilon_leaf', paramset_names)));
    a_root = paramset(find(strcmp('a_root', paramset_names)));
    a_stem = paramset(find(strcmp('a_stem', paramset_names)));
    b_root = paramset(find(strcmp('b_root', paramset_names)));
    b_stem = paramset(find(strcmp('b_stem', paramset_names)));
    
    y_theta_1 = y(1)/size_reservoir_stem_moles;
    y_theta_2 = y(2)/size_reservoir_leaf_moles;

    p_base = theta_to_p(y(1)/size_reservoir_stem_moles,...
        pi_0_stem, theta_r_stem, theta_tlp_stem, theta_ft_stem, epsilon_stem, S_stem); %, h_base, rhog_MPa);    
    
    p_leaf = theta_to_p(y(2)/size_reservoir_leaf_moles,...
        pi_0_leaf, theta_r_leaf, theta_tlp_leaf, theta_ft_leaf, epsilon_leaf, S_leaf); %, h_base, rhog_MPa);  
      
    flow_in_base_approx = K_max_root_total_moles*vc_integral_approx(p_base, p_soil, a_root, b_root)*...
        (p_soil - p_base - rhog_MPa* h_root/2) / (p_soil - p_base); % make sense to have h_root here be 0? no make average depth. Where is our reference?
    flow_in_base = vc_integral(p_base, p_soil, h_root/2, flow_in_base_approx,...
        K_max_root_total_moles, rhog_MPa, a_root, b_root);
    
    flow_out_base_approx = K_max_stem_moles*vc_integral_approx(p_leaf, p_base, a_stem, b_stem)*...
        (p_base - p_leaf - rhog_MPa* h_stem/2) / (p_base - p_leaf); % make sense to have h_root here be 0? no make average depth. Where is our reference?
    flow_out_base = vc_integral(p_leaf, p_base, h_stem/2, flow_out_base_approx,...
         K_max_stem_moles, rhog_MPa, a_stem, b_stem);   
    if t < 5 
        T = T_ini;
    elseif t >= 5 && t < 10
        T = 0.1*(t-5)+T_ini;
    else t >= 10
        T = 0.1*5+T_ini;
    end
    dydt = zeros(2,1);
    % dydt = zeros(1,1);
    dydt(1) = flow_in_base - flow_out_base;
    dydt(2) = flow_out_base - T;
    t
end
 
function [flux_approx] = vc_integral_approx(p_dos, p_ups, a, b)
%     if p_ups > 0 && p_dos > 0
%         error('positive pressure value in approx integral')
%     end
    t_dos = log(a .* exp(b.*p_dos) + 1);
    t_ups = log(a .* exp(b.*p_ups) + 1);
    flux_approx = (a+1) / (a*b) .* (t_ups - t_dos);
end

function [flux] = vc_integral(p_dos, p_ups, h, flux_approx, Kmax, rhog_MPa, a, b)
%     if p_ups <= 0 && p_dos <= 0
%         error('positive pressure value in approx integral')
%     end
    krghe = Kmax * rhog_MPa * h * (a+1) / a + flux_approx;
    lower = b * krghe;
    multi = Kmax * (a+1) / a .* flux_approx;
    upper_dos = log(a * krghe .* exp(b.*p_dos) + flux_approx);
    upper_ups = log(a * krghe .* exp(b.*p_ups) + flux_approx);
    flux = multi .* (upper_ups - upper_dos) ./ lower;
end
 
function p = theta_to_p(theta, pi_0, theta_r, theta_tlp, theta_ft, epsilon, S_s)      
    p = (theta-1)*5;  
end

function [theta] = p_to_theta(p, pi_0, pi_tlp, theta_r, theta_ft, epsilon, S_s)
    theta = p/5 + 1; 
end

% function p = theta_to_p(theta, pi_0, theta_r, theta_tlp, theta_ft, epsilon, S_s) 
%     if theta > theta_ft      
%         p = ((theta-theta_ft)./S_s)^2;
%     elseif theta > theta_tlp && theta <= theta_ft
%         p_sol = -abs(pi_0).*(theta_ft-theta_r)./(theta-theta_r);
%         p_p = abs(pi_0)-epsilon.*(theta_ft-theta)./(theta_ft-theta_r);
%         p = p_sol + p_p; 
%     elseif theta >= theta_r && theta <= theta_tlp
%         p_sol = -abs(pi_0).*(theta_ft-theta_r)./(theta-theta_r);
%         p = p_sol;
%     else
%         error('theta out of bounds')
%     end   
% end
% 
% function [theta] = p_to_theta(p, pi_0, pi_tlp, theta_r, theta_ft, epsilon, S_s)
%     if p > 0
%         theta = sqrt(p).*S_s+theta_ft;
%     elseif p > pi_tlp && p <= 0
%         theta = (epsilon*theta_r + epsilon*theta_ft - p.*theta_r + p.*theta_ft + theta_r*abs(pi_0) - theta_ft*abs(pi_0) - theta_r*(abs(pi_0).^2 + 2*epsilon.*p + 2*epsilon*abs(pi_0) - 2.*p.*abs(pi_0) + epsilon^2 + p.^2).^(1/2) + theta_ft*(abs(pi_0)^2 + 2*epsilon.*p + 2*epsilon*abs(pi_0) - 2.*p*abs(pi_0) + epsilon^2 + p.^2).^(1/2))/(2*epsilon);
%     elseif p <= pi_tlp
%         theta = theta_r + (abs(pi_0)*(theta_r - theta_ft))./p;
%     else 
%         error('p out of bounds')
%     end   
% end
 
% % Yujie's functions for v to p
% function p=theta_to_p(theta, theta_tlp, c_all, gas_R, T, theta_apo)
%     if theta > theta_tlp
%         p = -c_all * gas_R * T / (theta - theta_apo) + ...
%         epsilon_bulk * (theta - theta_tlp);
%     elseif theta > theta_apo
%         p = -c_all * gas_R * T / (theta - theta_apo);
%     else
%         error('theta out of bounds');
%     end
% end
% 
% function theta=p_to_theta(theta, theta_tlp, c_all, gas_R, T, theta_apo)
%     if p > p_tlp
%         theta = 
%         p = -c_all * gas_R * T / (theta - theta_apo) + ...
%         epsilon_bulk * (theta - theta_tlp);
%     elseif theta > theta_apo
%         
%         p = -c_all * gas_R * T / (theta - theta_apo);
%     else
%         error('theta out of bounds');
%     end
% end

