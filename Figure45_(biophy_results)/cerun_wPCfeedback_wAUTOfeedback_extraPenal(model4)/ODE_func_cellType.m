function dydt = ODE_func_cellType(t,y,RNAPHandle,p_cebp_hill,p_auto_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,remRate_constant_ratio_all_CEBPA,runx1_cellType_conc)

D_cebpa = log(2)/(100/60); % unit is 1/hour, same as PU1!!
tau_Mac = 40;

% pu1
dpudt = RNAPHandle(D7_valid_sublist,energy_table,[1,y(1),y(2),runx1_cellType_conc,1],nuc_destab_location)/...
    RNAPHandle(D7_valid_sublist,energy_table,[1,1,1,1,1],nuc_destab_location) - ... % production rate so far
    norm_remv_rate_const*y(1); % and removal rate

% cebpa
dcedt = cebpa_RegBy_pu1_Hill(p_cebp_hill,y(1))*cebpa_RegBy_cebpa_Hill(p_auto_hill,y(2))/(cebpa_RegBy_pu1_Hill(p_cebp_hill,1)*cebpa_RegBy_cebpa_Hill(p_auto_hill,1)) ...
    - remRate_constant_ratio_all_CEBPA*y(2); % production rate so far
%     ((log(2)/tau_cellType+D_cebpa)/(log(2)/tau_Mac+D_cebpa))*y(2);

dydt = [dpudt;dcedt];

