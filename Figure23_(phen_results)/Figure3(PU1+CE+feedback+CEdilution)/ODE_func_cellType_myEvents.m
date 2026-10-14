function [value,isterminal,direction] = ODE_func_cellType_myEvents(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,norm_remv_rate_const)

basal = pars_f1_hill(1);
vmax = pars_f1_hill(2) ;
N = pars_f1_hill(3);
khalf = pars_f1_hill(4);

basal_2 = pars_f2_hill(1);
vmax_2 = pars_f1_hill(2) ;
N_2 = pars_f1_hill(3);
khalf_2 = pars_f1_hill(4);

g_basal= pars_g_hill(1);
g_vmax = pars_g_hill(2);
g_N = pars_g_hill(3);
g_khalf = pars_g_hill(4);

D_pu1 = 1/70;
tau_Mac = 40;
% norm_remv_rate_const = (log(2)/tau_cellType+D_pu1)/(log(2)/tau_Mac+D_pu1);
norm_remv_rate_const_X = (log(2)/tau_cellType+D_X)/(log(2)/tau_Mac+D_X);

% function f (pu1)
dfdt = (basal+(vmax-basal)*(y(1)^N/(y(1)^N+khalf^N)))*...
    (basal_2+(vmax_2-basal_2)*(y(2)^N_2/(y(2)^N_2+khalf_2^N_2)))-...
    norm_remv_rate_const*y(1); 

% function g (X)
dgdt = g_basal+(g_vmax-g_basal)*(y(1)^g_N/(y(1)^g_N+g_khalf^g_N))-norm_remv_rate_const_X*y(2);

% construct the Jacobian matrix
% slope of a Hill function in terms of x: (v-b)*k^n*n*x^(n-1)/(k^n+x^n)^2
delfdelpu1 = (basal_2+(vmax_2-basal_2)*(y(2)^N_2/(y(2)^N_2+khalf_2^N_2)))*(vmax-basal)*khalf^N*N*y(1)^(N-1)/(khalf^N+y(1)^N)^2 - norm_remv_rate_const;
delfdelX = (basal+(vmax-basal)*(y(1)^N/(y(1)^N+khalf^N)))*(vmax_2-basal_2)*khalf_2^N_2*N_2*y(2)^(N_2-1)/(khalf_2^N_2+y(2)^N_2)^2;
delgdelpu1 = (g_vmax-g_basal)*g_khalf^g_N*g_N*y(1)^(g_N-1)/(g_khalf^g_N+y(1)^g_N)^2;
delgdelX = -norm_remv_rate_const_X;

J = [delfdelpu1,delfdelX;delgdelpu1,delgdelX];
eig_ODE = eig(J);
bool_sign_real_eig = all(real(eig_ODE)<0); % real(eig_ODE)<0 will give a 1 
% or 0 array indicate whehter each eigenvalues is < 0. if product of all elements in this array is 1, then return true
bool_equilibrium = (abs(dgdt) < 1e-5)*(abs(dfdt) < 1e-5);
% if the change in dgdt and dfdt are small enough

value = bool_equilibrium*bool_sign_real_eig - 1; % if value = 0, stop integration
isterminal = 1;
direction = 0;
end