function dydt = ODE_func_cellType(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,tau_cellType,norm_remv_rate_const_PU1)

basal = pars_f1_hill(1);
vmax = pars_f1_hill(2) ;
N = pars_f1_hill(3);
khalf = pars_f1_hill(4);

basal_2 = pars_f2_hill(1);
vmax_2 = pars_f2_hill(2) ;
N_2 = pars_f2_hill(3);
khalf_2 = pars_f2_hill(4);

g_basal= pars_g_hill(1);
g_vmax = pars_g_hill(2);
g_N = pars_g_hill(3);
g_khalf = pars_g_hill(4);

% D_pu1 = 1/70;
tau_Mac = 40;
% norm_remv_rate_const_PU1 = (log(2)/tau_cellType+D_pu1)/(log(2)/tau_Mac+D_pu1);
norm_remv_rate_const_X = (log(2)/tau_cellType+D_X)/(log(2)/tau_Mac+D_X);

% function f (pu1)
dfdt = (basal+(vmax-basal)*(y(1)^N/(y(1)^N+khalf^N)))*...
    (basal_2+(vmax_2-basal_2)*(y(2)^N_2/(y(2)^N_2+khalf_2^N_2)))-...
    norm_remv_rate_const_PU1*y(1); 

% function g (X)
dgdt = g_basal+(g_vmax-g_basal)*(y(1)^g_N/(y(1)^g_N+g_khalf^g_N))-norm_remv_rate_const_X*y(2);

dydt = [dfdt;dgdt];

