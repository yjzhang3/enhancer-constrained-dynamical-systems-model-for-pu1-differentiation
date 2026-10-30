function [err,tB,tPL,tPH,tmac,...
    Bcell_over_time,PL_over_time,PH_over_time,Mac_over_time,...
    pred_pu1_ss,pred_pu1_prod] = obj_twoHill_ODE(p)

load(fullfile(fileparts(pwd), 'Error propagation and Normalization', 'error_prop.mat'));

p_Hill = p(1:8);
basal = p_Hill(1);
vmax = p_Hill(2) ;
N = p_Hill(3);
khalf = p_Hill(4);
g_basal = p_Hill(5);
g_vmax = p_Hill(6);
g_N = p_Hill(7);
g_khalf = p_Hill(8);
g_input_Bcell = p(9);
g_input_ProLow = p(10);
g_input_ProHigh = p(11);
g_input_mac = p(12);

% function setup
ODE_func_Bcell = @(t,y) ODE_func(t,y,p_Hill,g_input_Bcell,remRate_constant_ratio_all(1,2));
ODE_func_ProLow = @(t,y) ODE_func(t,y,p_Hill,g_input_ProLow,remRate_constant_ratio_all(2,2));
ODE_func_ProHigh = @(t,y) ODE_func(t,y,p_Hill,g_input_ProHigh,remRate_constant_ratio_all(3,1));
ODE_func_Mac = @(t,y) ODE_func(t,y,p_Hill,g_input_mac,1);

tspan = [0 50];
[tB,Bcell_over_time] = ode45(ODE_func_Bcell,tspan,[conc_ratio_all(1,2)],odeset(RelTol=1e-5));
[tPL,PL_over_time] = ode45(ODE_func_ProLow,tspan,[conc_ratio_all(1,2)],odeset(RelTol=1e-5));
[tPH,PH_over_time] = ode45(ODE_func_ProHigh,tspan,[conc_ratio_all(1,2)],odeset(RelTol=1e-5));
[tmac,Mac_over_time] = ode45(ODE_func_Mac,tspan,[conc_ratio_all(1,2)],odeset(RelTol=1e-5));

% predict pu1 ss in diff ct
pred_pu1_ss = [Bcell_over_time(end,1),PL_over_time(end,1),PH_over_time(end,1),Mac_over_time(end,1)];
real_pu1_ss = [conc_ratio_all(1,2),conc_ratio_all(2,2),conc_ratio_all(3,2),1];
er1 = weighted_msd(pred_pu1_ss,real_pu1_ss);

% predict production rate
production_rate = @(x,w) (basal+(vmax-basal)*(x^N/(x^N+khalf^N)))*g_basal+(g_vmax-g_basal)*(w^g_N/(w^g_N+g_khalf^g_N));
pred_pu1_prod = [production_rate(pred_pu1_ss(1),g_input_Bcell),...
    production_rate(pred_pu1_ss(2),g_input_ProLow),...
    production_rate(pred_pu1_ss(3),g_input_ProHigh),...
    production_rate(pred_pu1_ss(4),g_input_mac)];
real_pu1_prod = [prodRate_ratio_all(1,2),prodRate_ratio_all(2,2),prodRate_ratio_all(3,2),1];
er2 = weighted_msd(pred_pu1_prod,real_pu1_prod);


% error (difference between prod and remv)
err = er1+er2;
end

function dydt = ODE_func(t,y,p_Hill,w_ct,norm_remv_rate_cellType)

basal = p_Hill(1);
vmax = p_Hill(2) ;
N = p_Hill(3);
khalf = p_Hill(4);
g_basal = p_Hill(5);
g_vmax = p_Hill(6);
g_N = p_Hill(7);
g_khalf = p_Hill(8);

production_rate = (basal+(vmax-basal)*(y(1)^N/(y(1)^N+khalf^N)))*(g_basal+(g_vmax-g_basal)*(w_ct^g_N/(w_ct^g_N+g_khalf^g_N)));
removal_rate = norm_remv_rate_cellType*y(1);
dydt = production_rate - removal_rate;

end

