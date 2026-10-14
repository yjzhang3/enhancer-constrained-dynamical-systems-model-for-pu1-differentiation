function [final_err,pred_pu1_ss,pred_X_ss,pred_pu1_prod,...
    tb,tpl,tph,tm,...
    bcell_over_time,prolow_over_time,prohigh_over_time,mac_over_time] = obj_twoHill_wSlope_wFeedback(p)

load('error_prop_0216.mat');

% b, v, N, K
pars_f1_hill = p(1:4); % autoregulation
pars_f2_hill = p(5:8); % PU1 regulated by X
pars_g_hill = p(9:12); % X regulated by PU1
D_X = p(13); % natural degradation rate of TF X

prodHandle_raw = @(pars_hill,x) pars_hill(1)+(pars_hill(2)-pars_hill(1))*(x^pars_hill(3)/(x^pars_hill(3)+pars_hill(4)^pars_hill(3)));
prodHandle_final = @(x1,x2) prodHandle_raw(pars_f1_hill,x1)*prodHandle_raw(pars_f2_hill,x2);

tau_Mac = [37.8,40,56.4];
tau_Bcell = [17.7,25.5,41.2];
tau_ProLow = [19.9,28.8,47.8];
tau_ProHigh = [9.47,11.6,15.8];

%% function setup (ct specific) and solve for ss
tspan = [0 400];

% tau_cellType is to computedilution rate for X inly
bcell_ode_fun = @(t,y) ODE_func_cellType(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,tau_Bcell(2),remRate_constant_ratio_all(1,2));
prolow_ode_fun = @(t,y) ODE_func_cellType(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,tau_ProLow(2),remRate_constant_ratio_all(2,2));
prohigh_ode_fun = @(t,y) ODE_func_cellType(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,tau_ProHigh(2),remRate_constant_ratio_all(3,1));
mac_ode_fun = @(t,y) ODE_func_cellType(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,tau_Mac(2),1);


%%% with events ode45(func,tspan,y0,options)
% bcell_event = @(t,y) ODE_func_cellType_myEvents(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,remRate_constant_ratio_all(1,2));
% prolow_event = @(t,y) ODE_func_cellType_myEvents(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,remRate_constant_ratio_all(2,2));
% prohigh_event = @(t,y) ODE_func_cellType_myEvents(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,remRate_constant_ratio_all(3,1));
% mac_event = @(t,y) ODE_func_cellType_myEvents(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,1);
% [tb,bcell_over_time,te,ye,~] = ode45(bcell_ode_fun,tspan,[conc_ratio_all(1,2),conc_ratio_all(1,2)],odeset('Events',bcell_event,'RelTol',1e-6));
% [tpl,prolow_over_time,te,ye,~] = ode45(prolow_ode_fun,tspan,[conc_ratio_all(2,2),conc_ratio_all(2,2)],odeset('Events',prolow_event,'RelTol',1e-6));
% [tph,prohigh_over_time,te,ye,~] = ode45(prohigh_ode_fun,tspan,[conc_ratio_all(3,2),conc_ratio_all(3,2)],odeset('Events',prohigh_event,'RelTol',1e-6));
% [tm,mac_over_time,te,ye,~] = ode45(mac_ode_fun,tspan,[1,1],odeset('Events',mac_event));


%%% without events
[tb,bcell_over_time] = ode45(bcell_ode_fun,tspan,[conc_ratio_all(1,2),conc_ratio_all(1,2)],odeset('RelTol',1e-5));
[tpl,prolow_over_time] = ode45(prolow_ode_fun,tspan,[conc_ratio_all(2,2),conc_ratio_all(2,2)],odeset('RelTol',1e-5));
[tph,prohigh_over_time] = ode45(prohigh_ode_fun,tspan,[conc_ratio_all(3,2),conc_ratio_all(3,2)],odeset('RelTol',1e-5));
[tm,mac_over_time] = ode45(mac_ode_fun,tspan,[1,1],odeset('RelTol',1e-5));

pred_pu1_ss = [bcell_over_time(end,1),prolow_over_time(end,1),prohigh_over_time(end,1),mac_over_time(end,1)];
pred_X_ss = [bcell_over_time(end,2),prolow_over_time(end,2),prohigh_over_time(end,2),mac_over_time(end,2)];
pred_pu1_prod = [prodHandle_final(pred_pu1_ss(1),pred_X_ss(1)),...
    prodHandle_final(pred_pu1_ss(2),pred_X_ss(2)),...
    prodHandle_final(pred_pu1_ss(3),pred_X_ss(3)),...
    prodHandle_final(pred_pu1_ss(4),pred_X_ss(4))];

%% enforce steady state of PU1 of each cell type happens at data-constrained range
% Bcell
g1 = 0.1;
g2 = 0.1;
ksmall = conc_ratio_all(1,1);
klarge = pred_pu1_ss(1);
p1 = penalty_function(ksmall,klarge,g1,g2);

ksmall = pred_pu1_ss(1);
klarge = conc_ratio_all(1,3);
p2 = penalty_function(ksmall,klarge,g1,g2);

% Prolow
g1 = 0.1;
g2 = 0.1;
ksmall = conc_ratio_all(2,1);
klarge = pred_pu1_ss(2);
p3 = penalty_function(ksmall,klarge,g1,g2);

ksmall = pred_pu1_ss(2);
klarge = conc_ratio_all(2,3);
p4 = penalty_function(ksmall,klarge,g1,g2);

% ProHigh
g1 = 0.1;
g2 = 0.1;
ksmall = conc_ratio_all(3,1);
klarge = pred_pu1_ss(3);
p5 = penalty_function(ksmall,klarge,g1,g2);

ksmall = pred_pu1_ss(3);
klarge = conc_ratio_all(3,3);
p6 = penalty_function(ksmall,klarge,g1,g2);

%% enforce production lies in data range
% Bcell
g1 = 0.1;
g2 = 0.1;
ksmall = prodRate_ratio_all(1,1);
klarge = pred_pu1_prod(1);
p21 = penalty_function(ksmall,klarge,g1,g2);

ksmall = pred_pu1_prod(1);
klarge = prodRate_ratio_all(1,3);
p22 = penalty_function(ksmall,klarge,g1,g2);

% Prolow
g1 = 0.1;
g2 = 0.1;
ksmall = prodRate_ratio_all(2,1);
klarge = pred_pu1_prod(2);
p23 = penalty_function(ksmall,klarge,g1,g2);

ksmall = pred_pu1_prod(2);
klarge = prodRate_ratio_all(2,3);
p24 = penalty_function(ksmall,klarge,g1,g2);

% ProHigh
g1 = 0.1;
g2 = 0.1;
ksmall = prodRate_ratio_all(3,1);
klarge = pred_pu1_prod(3);
p25 = penalty_function(ksmall,klarge,g1,g2);

% ksmall = pred_pu1_prod(3);
% klarge = prodRate_ratio_all(3,3);
% p26 = penalty_function(ksmall,klarge,g1,g2);

%% macrophage is always 1
err = weighted_msd([pred_pu1_ss(4),pred_pu1_prod(4)],[1,1]);

final_err = p1 + p2 + p3 + p4 + p5 + p6 + p21 + p22 + p23 + p24 + p25 + 0 + err;

end

