function [tb,Bcell_conc_over_time,tpl,ProLow_conc_over_time,tph,ProHigh_conc_over_time] = ...
    debug_ODE_stability(best_pu1_ss,best_cebpa_ss,indep_predicted_p_cebp_hill,energy_table,runx1_ss_conc)
%% housekeeping
load('error_prop_0216.mat'); 
% load('ct_TFPROD_constraint.mat',"ce_Bcell_over_mac_bound","ce_pro_over_mac_bound")

RNAPHandle = @RNAP_occ_prob_cebpa_runx1_general;

selec = 1:12;
nuc_destab_location = [0,1,0,1,0,0];
site_index = 1;

runx1_bcell_conc=runx1_ss_conc(1);
runx1_ProLow_conc=runx1_ss_conc(2);
runx1_ProHigh_conc=runx1_ss_conc(3);

% D7 
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D7_valid_sublist = valid_sublist;

%% construct cell-type specific ODE
% remember D7 is wildtype whose production we're interested in computing
% concentration table order of entry: [nuc],[pu1],[cebp],[runx],[RNAP]

ODE_func_Bcell = @(t,y) debug_ODE_func_cellType_noAUTO(t,y,RNAPHandle,indep_predicted_p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(1,2),remRate_constant_ratio_all_CEBPA(1,2),runx1_bcell_conc);
ODE_func_ProLow = @(t,y) debug_ODE_func_cellType_noAUTO(t,y,RNAPHandle,indep_predicted_p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(2,2),remRate_constant_ratio_all_CEBPA(2,2),runx1_ProLow_conc);
ODE_func_ProHigh = @(t,y) debug_ODE_func_cellType_noAUTO(t,y,RNAPHandle,indep_predicted_p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(3,1),remRate_constant_ratio_all_CEBPA(3,1),runx1_ProHigh_conc);

% fs_func_Bcell = @(y) ODE_func_cellType([],y,RNAPHandle,p_cebp_hill,p_auto_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(1,2),remRate_constant_ratio_all_CEBPA(1,2),runx1_bcell_conc);
% fs_func_ProLow = @(y) ODE_func_cellType([],y,RNAPHandle,p_cebp_hill,p_auto_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(2,2),remRate_constant_ratio_all_CEBPA(2,2),runx1_ProLow_conc);
% fs_func_ProHigh = @(y) ODE_func_cellType([],y,RNAPHandle,p_cebp_hill,p_auto_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(3,1),remRate_constant_ratio_all_CEBPA(3,1),runx1_ProHigh_conc);

%% quickly compute an equilibrium point using the ODE solver
tspan = 200;
[tb,Bcell_conc_over_time] = ode23t(ODE_func_Bcell,[0 tspan],[best_pu1_ss(1),best_cebpa_ss(1)]',odeset('RelTol',1e-5,'AbsTol',1e-6));
[tpl,ProLow_conc_over_time] = ode23t(ODE_func_ProLow,[0 tspan],[best_pu1_ss(2),best_cebpa_ss(2)]',odeset('RelTol',1e-5,'AbsTol',1e-6));
[tph,ProHigh_conc_over_time] = ode23t(ODE_func_ProHigh,[0 tspan],[best_pu1_ss(3),best_cebpa_ss(3)]',odeset('RelTol',1e-5,'AbsTol',1e-6));



