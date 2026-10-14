function [dydt] = ODE_func_cellType_wCEAUTO(t,y,RNAPHandle,p_cebp_hill,p_auto_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,norm_remv_rate_const_CEBPA)

% pu1
dpudt = RNAPHandle(D7_valid_sublist,energy_table,[1,y(1),y(2),1],nuc_destab_location)/...
    RNAPHandle(D7_valid_sublist,energy_table,[1,1,1,1],nuc_destab_location) - ... % production rate so far
    norm_remv_rate_const*y(1); % and removal rate

% pu1_prod = RNAPHandle(D7_valid_sublist,energy_table,[1,y(1),y(2),1],nuc_destab_location)/...
%     RNAPHandle(D7_valid_sublist,energy_table,[1,1,1,1],nuc_destab_location);
% pu1_remv = norm_remv_rate_const*y(1);

% cebpa
dcedt = cebpa_RegBy_pu1_Hill(p_cebp_hill,y(1))*cebpa_RegBy_cebpa_Hill(p_auto_hill,y(2))/(cebpa_RegBy_pu1_Hill(p_cebp_hill,1)*cebpa_RegBy_cebpa_Hill(p_auto_hill,1)) ...
    - norm_remv_rate_const_CEBPA*y(2);
    % ((log(2)/tau_cellType+D_cebpa)/(log(2)/tau_Mac+D_cebpa))*y(2);
    
dydt = [dpudt;dcedt];

