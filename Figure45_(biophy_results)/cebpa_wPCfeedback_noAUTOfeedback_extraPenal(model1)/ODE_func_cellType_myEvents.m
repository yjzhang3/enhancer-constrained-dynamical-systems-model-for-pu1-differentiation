function [value,isterminal,direction] = ODE_func_cellType_myEvents(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,tau_cellType,norm_remv_rate_const,norm_remv_rate_const_CEBPA)

% pu1
dpudt = RNAPHandle(D7_valid_sublist,energy_table,[1,y(1),y(2),1],nuc_destab_location)/...
    RNAPHandle(D7_valid_sublist,energy_table,[1,1,1,1],nuc_destab_location) - ... % production rate so far
    norm_remv_rate_const*y(1); % and removal rate

% pu1_prod = RNAPHandle(D7_valid_sublist,energy_table,[1,y(1),y(2),1],nuc_destab_location)/...
%     RNAPHandle(D7_valid_sublist,energy_table,[1,1,1,1],nuc_destab_location);
% pu1_remv = norm_remv_rate_const*y(1);

% cebpa
dcedt = cebpa_RegBy_pu1_Hill(p_cebp_hill,y(1))/cebpa_RegBy_pu1_Hill(p_cebp_hill,1) - norm_remv_rate_const_CEBPA*y(2);
    % ((log(2)/tau_cellType+D_cebpa)/(log(2)/tau_Mac+D_cebpa))*y(2);
% construct the Jacobian matrix

% slope of a Hill function in terms of x: (v-b)*k^n*n*x^(n-1)/(k^n+x^n)^2
% delfdelpu1 = (basal_2+(vmax_2-basal_2)*(y(2)^N_2/(y(2)^N_2+khalf_2^N_2)))*(vmax-basal)*khalf^N*N*y(1)^(N-1)/(khalf^N+y(1)^N)^2 - norm_remv_rate_const;
% delfdelX = (basal+(vmax-basal)*(y(1)^N/(y(1)^N+khalf^N)))*(vmax_2-basal_2)*khalf_2^N_2*N_2*y(2)^(N_2-1)/(khalf_2^N_2+y(2)^N_2)^2;
% delgdelpu1 = (g_vmax-g_basal)*g_khalf^g_N*g_N*y(1)^(g_N-1)/(g_khalf^g_N+y(1)^g_N)^2;
% delgdelX = -norm_remv_rate_const_X;

% J = [delfdelpu1,delfdelX;delgdelpu1,delgdelX];
% eig_ODE = eig(J);
% bool_sign_real_eig = all(real(eig_ODE)<0); % real(eig_ODE)<0 will give a 1 
% or 0 array indicate whehter each eigenvalues is < 0. if product of all elements in this array is 1, then return true
bool_equilibrium = (abs(dpudt) < 1e-8)*(abs(dcedt) < 1e-8);
% if the change in dgdt and dfdt are small enough

% value = bool_equilibrium*bool_sign_real_eig - 1; % if value = 0, stop integration
value = bool_equilibrium - 1;
isterminal = 1;
direction = 0;

end