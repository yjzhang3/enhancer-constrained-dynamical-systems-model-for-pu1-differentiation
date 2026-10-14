function [final_err,p_cebp_hill,recomputed_concentration_table_Bcell,recomputed_concentration_table_ProLow,recomputed_concentration_table_ProHigh,...
    real_d_final,sim_d_final,D7_valid_sublist] = ...
    ODE_objective_debug_model3(selec,p,nuc_destab_location,energy_table,promoter_only_rescale_factor,pu1_ss_conc,cebpa_ss_conc,runx1_ss_conc)
% directly modified from objective_with_cebpa_general_PU1CEfeedback, so see
% information in the other file

%% housekeeping
load('error_prop_0216.mat'); 
load('ct_TFPROD_constraint.mat',"ce_Bcell_over_mac_bound","ce_pro_over_mac_bound")

RNAPHandle = @RNAP_occ_prob_cebpa_runx1_general;
site_index = 1;

% D7 
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D7_valid_sublist = valid_sublist;

%% load pars
% G_NUC = p(1);
% G_PU1_14 = p(2);
% G_PU1_12 = p(3);
% G_PU1_10 = p(4);
% G_CEBP = p(5);
% G_RUNX = p(6);
% G_RNAP = p(7);
% G_NUC_PU1 = p(8); 
% G_NUC_CEBP = p(9);
% G_NUC_RUNX1 = p(10);
% G_PU1_RNAP = p(11);
% G_CEBP_RNAP = p(12);
% G_RUNX1_RNAP = p(13);
% G_RUNX_PU1 = p(14);
% G_RUNX_CEBP = p(15);
% G_PU1_CEBP = p(16);

% promoter_only_rescale_factor = p(17);
runx1_bcell_conc = runx1_ss_conc(1);
runx1_ProLow_conc = runx1_ss_conc(2);
runx1_ProHigh_conc = runx1_ss_conc(3);

p_cebp_hill = p(1:4);

pu1_Bcell_ss_conc = pu1_ss_conc(1);
pu1_ProLow_ss_conc = pu1_ss_conc(2);
pu1_ProHigh_ss_conc = pu1_ss_conc(3);

cebpa_bcell_ss_conc = cebpa_ss_conc(1);
cebpa_ProLow_ss_conc = cebpa_ss_conc(2);
cebpa_ProHigh_ss_conc = cebpa_ss_conc(3);

% energy_table = [G_NUC,G_PU1_14,G_PU1_12,G_PU1_10,G_CEBP,G_RUNX,G_RNAP, ... %1-7
%     G_NUC_PU1,G_NUC_CEBP,G_NUC_RUNX1,... % 8-10
%     G_PU1_RNAP,G_CEBP_RNAP,G_RUNX1_RNAP,... % 11-13
%     G_RUNX_PU1,G_RUNX_CEBP,G_PU1_CEBP]; ... % 14-16

% nuc_destab_location =
% [PU1_PRO,PU1_allTFBS,C/EBPa_PRO,C/EBPa_12,RUNX1_PRO,RUNX1_12]

%% construct cell-type specific ODE
% remember D7 is wildtype whose production we're interested in computing
% concentration table order of entry: [nuc],[pu1],[cebp],[runx],[RNAP]

ODE_func_Bcell = @(t,y) ODE_func_cellType([],y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(1,2),remRate_constant_ratio_all_CEBPA(1,2),runx1_bcell_conc);
ODE_func_ProLow = @(t,y) ODE_func_cellType([],y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(2,2),remRate_constant_ratio_all_CEBPA(2,2),runx1_ProLow_conc);
ODE_func_ProHigh = @(t,y) ODE_func_cellType([],y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(3,1),remRate_constant_ratio_all_CEBPA(3,1),runx1_ProHigh_conc);

fs_func_Bcell = @(y) ODE_func_cellType([],y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(1,2),remRate_constant_ratio_all_CEBPA(1,2),runx1_bcell_conc);
fs_func_ProLow = @(y) ODE_func_cellType([],y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(2,2),remRate_constant_ratio_all_CEBPA(2,2),runx1_ProLow_conc);
fs_func_ProHigh = @(y) ODE_func_cellType([],y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(3,1),remRate_constant_ratio_all_CEBPA(3,1),runx1_ProHigh_conc);
% 
%% quickly compute an equilibrium point using the ODE solver
tspan = 200;
[~,Bcell_conc_over_time] = ode23t(ODE_func_Bcell,[0 tspan],[pu1_Bcell_ss_conc,cebpa_bcell_ss_conc]',odeset('RelTol',1e-5,'AbsTol',1e-6));
[~,ProLow_conc_over_time] = ode23t(ODE_func_ProLow,[0 tspan],[pu1_ProLow_ss_conc,cebpa_ProLow_ss_conc]',odeset('RelTol',1e-5,'AbsTol',1e-6));
[~,ProHigh_conc_over_time] = ode23t(ODE_func_ProHigh,[0 tspan],[pu1_ProHigh_ss_conc,cebpa_ProHigh_ss_conc]',odeset('RelTol',1e-5,'AbsTol',1e-6));

%% then use fsolve (implicit solver)
options = optimset('Display','off');
Bcell_sol = fsolve(fs_func_Bcell,Bcell_conc_over_time(end,:),options);
ProLow_sol = fsolve(fs_func_ProLow,ProLow_conc_over_time(end,:),options);
ProHigh_sol = fsolve(fs_func_ProHigh,ProHigh_conc_over_time(end,:),options);

%% test stability using eigenvalues
Bcell_sol_stab = zeros(length(Bcell_sol(:,1)),2);
for i = 1:length(Bcell_sol(:,1))
    Bcell_sol_stab(i,:) = compute_Jacobian_eig_numerical(ODE_func_Bcell,[Bcell_sol(i,1),Bcell_sol(i,2)]);
end
idx1 = find(max(Bcell_sol_stab,[],2) < 0); % find stable sols

ProLow_sol_stab = zeros(length(ProLow_sol(:,1)),2);
for i = 1:length(ProLow_sol(:,1))
    ProLow_sol_stab(i,:) = compute_Jacobian_eig_numerical(ODE_func_ProLow,[ProLow_sol(i,1),ProLow_sol(i,2)]);
end
idx2 = find(max(ProLow_sol_stab,[],2) < 0);

ProHigh_sol_stab = zeros(length(ProHigh_sol(:,1)),2);
for i = 1:length(ProHigh_sol(:,1))
    ProHigh_sol_stab(i,:) = compute_Jacobian_eig_numerical(ODE_func_ProHigh,[ProHigh_sol(i,1),ProHigh_sol(i,2)]);
end
idx3 = find(max(ProHigh_sol_stab,[],2) < 0);

%% extract steady state
if ~isempty(idx1)
    recomputed_pu1_Bcell_ss_conc = Bcell_sol(idx1,1);
    recomputed_cebpa_bcell_ss_conc = Bcell_sol(idx1,2);
else
    [tB2,Bcell_over_time2] = ode23t(ODE_func_Bcell,[0 tspan*2],Bcell_sol(end,:),odeset('RelTol',1e-10,'AbsTol',1e-12));
    recomputed_pu1_Bcell_ss_conc = Bcell_over_time2(end,1);
    recomputed_cebpa_bcell_ss_conc = Bcell_over_time2(end,2);
end

if ~isempty(idx2)
    recomputed_pu1_ProLow_ss_conc = ProLow_sol(idx2,1);
    recomputed_cebpa_ProLow_ss_conc = ProLow_sol(idx2,2);
else
    [tPL2,ProLow_over_time2] = ode23t(ODE_func_Bcell,[0 tspan*2],ProLow_sol(end,:),odeset('RelTol',1e-10,'AbsTol',1e-12));
    recomputed_pu1_ProLow_ss_conc = ProLow_over_time2(end,1);
    recomputed_cebpa_ProLow_ss_conc = ProLow_over_time2(end,2);
end

if ~isempty(idx3)
    recomputed_pu1_ProHigh_ss_conc = ProHigh_sol(idx3,1);
    recomputed_cebpa_ProHigh_ss_conc = ProHigh_sol(idx3,2);
else
    [tPH2,ProHigh_over_time2] = ode23t(ODE_func_ProHigh,[0 tspan*2],ProHigh_sol(end,:),odeset('RelTol',1e-10,'AbsTol',1e-12));
    recomputed_pu1_ProHigh_ss_conc = ProHigh_over_time2(end,1);
    recomputed_cebpa_ProHigh_ss_conc = ProHigh_over_time2(end,2);
end
% 
recomputed_pu1_ss_conc = [recomputed_pu1_Bcell_ss_conc,recomputed_pu1_ProLow_ss_conc,recomputed_pu1_ProHigh_ss_conc,1];

% cebpa_ss_conc = [cebpa_bcell_ss_conc,cebpa_ProLow_ss_conc,cebpa_ProHigh_ss_conc,1];

% runx1_ss_conc = [runx1_bcell_conc,runx1_ProLow_conc,runx1_ProHigh_conc,1];

concentration_table_Mac = [1,1,1,1,1]; % [nuc],[pu1],[cebp],[RNAP] for *macrophage*
% concentration_table_Bcell = [1,pu1_Bcell_ss_conc,cebpa_bcell_ss_conc,runx1_bcell_conc,1]; % [nuc],[pu1],[cebp],[RNAP]
% concentration_table_ProLow = [1,pu1_ProLow_ss_conc,cebpa_ProLow_ss_conc,runx1_ProLow_conc,1];
% concentration_table_ProHigh = [1,pu1_ProHigh_ss_conc,cebpa_ProHigh_ss_conc,runx1_ProHigh_conc,1]; 

recomputed_concentration_table_Mac = concentration_table_Mac;
recomputed_concentration_table_Bcell = [1,recomputed_pu1_Bcell_ss_conc,recomputed_cebpa_bcell_ss_conc,runx1_bcell_conc,1]; % [nuc],[pu1],[cebp],[RNAP]
recomputed_concentration_table_ProLow = [1,recomputed_pu1_ProLow_ss_conc,recomputed_cebpa_ProLow_ss_conc,runx1_ProLow_conc,1];
recomputed_concentration_table_ProHigh = [1,recomputed_pu1_ProHigh_ss_conc,recomputed_cebpa_ProHigh_ss_conc,runx1_ProHigh_conc,1]; 

%% enforce steady state of PU1 of each cell type happens at data-constrained range
% Bcell
g1 = 1;
g2 = 0.1;
ksmall = conc_ratio_all(1,1);
klarge = recomputed_pu1_Bcell_ss_conc;
p1 = penalty_function(ksmall,klarge,g1,g2);

ksmall = recomputed_pu1_Bcell_ss_conc;
klarge = conc_ratio_all(1,3);
p2 = penalty_function(ksmall,klarge,g1,g2);

% Prolow
ksmall = conc_ratio_all(2,1);
klarge = recomputed_pu1_ProLow_ss_conc;
p3 = penalty_function(ksmall,klarge,g1,g2);

ksmall = recomputed_pu1_ProLow_ss_conc;
klarge = conc_ratio_all(2,3);
p4 = penalty_function(ksmall,klarge,g1,g2);

% ProHigh

ksmall = conc_ratio_all(3,1);
klarge = recomputed_pu1_ProHigh_ss_conc;
p5 = penalty_function(ksmall,klarge,g1,g2);

ksmall = recomputed_pu1_ProHigh_ss_conc;
klarge = conc_ratio_all(3,3);
p6 = penalty_function(ksmall,klarge,g1,g2);

% pu1ss_error = [p1,p2,p3,p4,p5,p6];

%% enforce CEBPa production rate lies in data-constrained range
% % Bcell
% 
% ksmall = ce_Bcell_over_mac_bound(1);
% klarge = cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_Bcell_ss_conc)/cebpa_RegBy_pu1_Hill(p_cebp_hill,1);
% p11 = penalty_function(ksmall,klarge,g1,g2);
% 
% ksmall = cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_Bcell_ss_conc)/cebpa_RegBy_pu1_Hill(p_cebp_hill,1);
% klarge = ce_Bcell_over_mac_bound(2);
% p12 = penalty_function(ksmall,klarge,g1,g2);
% 
% % Prolow
% 
% ksmall = ce_pro_over_mac_bound(1);
% klarge = cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_ProLow_ss_conc)/cebpa_RegBy_pu1_Hill(p_cebp_hill,1);
% p13 = penalty_function(ksmall,klarge,g1,g2);
% 
% ksmall = cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_ProLow_ss_conc)/cebpa_RegBy_pu1_Hill(p_cebp_hill,1);
% klarge = ce_pro_over_mac_bound(2);
% p14 = penalty_function(ksmall,klarge,g1,g2);
% 
% % ProHigh
% 
% ksmall = ce_pro_over_mac_bound(1);
% klarge = cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_ProHigh_ss_conc)/cebpa_RegBy_pu1_Hill(p_cebp_hill,1);
% p15 = penalty_function(ksmall,klarge,g1,g2);
% 
% ksmall = cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_ProHigh_ss_conc)/cebpa_RegBy_pu1_Hill(p_cebp_hill,1);
% klarge = ce_pro_over_mac_bound(2);
% p16 = penalty_function(ksmall,klarge,g1,g2);
% 
% cebpaprod_error = [p11,p12,p13,p14,p15,p16];

%% enforce production lies in data range
% Bcell
recomputed_D7_prob_bcell = RNAPHandle(D7_valid_sublist,energy_table,recomputed_concentration_table_Bcell,nuc_destab_location)/...
    RNAPHandle(D7_valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

ksmall = prodRate_ratio_all(1,1);
klarge = recomputed_D7_prob_bcell;
p21 = penalty_function(ksmall,klarge,g1,g2);

ksmall = recomputed_D7_prob_bcell;
klarge = prodRate_ratio_all(1,3);
p22 = penalty_function(ksmall,klarge,g1,g2);

% Prolow
recomputed_D7_prob_ProLow = RNAPHandle(D7_valid_sublist,energy_table,recomputed_concentration_table_ProLow,nuc_destab_location)/...
    RNAPHandle(D7_valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

ksmall = prodRate_ratio_all(2,1);
klarge = recomputed_D7_prob_ProLow;
p23 = penalty_function(ksmall,klarge,g1,g2);

ksmall = recomputed_D7_prob_ProLow;
klarge = prodRate_ratio_all(2,3);
p24 = penalty_function(ksmall,klarge,g1,g2);

% ProHigh
recomputed_D7_prob_ProHigh = RNAPHandle(D7_valid_sublist,energy_table,recomputed_concentration_table_ProHigh,nuc_destab_location)/...
    RNAPHandle(D7_valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

ksmall = prodRate_ratio_all(3,1);
klarge = recomputed_D7_prob_ProHigh;
p25 = penalty_function(ksmall,klarge,g1,g2);

ksmall = recomputed_D7_prob_ProHigh;
klarge = prodRate_ratio_all(3,3);
p26 = penalty_function(ksmall,klarge,g1,g2);

% pu1_ss_prod = [D7_prob_bcell,D7_prob_ProLow,D7_prob_ProHigh];
% pu1prod_error = [p21,p22,p23,p24,p25,p26];
%% START: make sure production rate with or without mutation is consistent with PU1 gene expression data
%%%%%%%%%%%%%%%%%%%%%%

%% real data
ratio_D2_over_D1 = 1.6; %1
ratio_D3_over_D1 = 7.06; %2
ratio_D4_over_D1 = 229; %3
ratio_D5_over_D1 = 275.4; %4
ratio_D6_over_D1 = 2556.3; %5
ratio_D7_over_D1 = 1471.3; %6
ratio_F2_over_D1 = 286.30; %7
ratio_F3_over_D1 = 237.73; %8
ratio_F4_over_D1 = 35.27; %9
ratio_E2_over_D1 = 26.7*promoter_only_rescale_factor; %10
ratio_E3_over_D1 = 20.3*promoter_only_rescale_factor; %11

ratio_AML1_paper = ratio_D4_over_D1*32/880; %12 % without AML1 binding sites, decrease in PU1 expression
% ratio_CEBP_paper_C1mut = 77/100; % do you really need all CEBPa binding sites? The effect of taking out all 3 cebap sites and taking out the whole -14kb region is the same!
% ratio_CEBP_paper_C2mut = 59/100;
% ratio_CEBP_paper_C5mut = 20/100;

real_d = [ratio_D2_over_D1,ratio_D3_over_D1,ratio_D4_over_D1,ratio_D5_over_D1,ratio_D6_over_D1,ratio_D7_over_D1,ratio_F2_over_D1,ratio_F3_over_D1,ratio_F4_over_D1,ratio_E2_over_D1,ratio_E3_over_D1,ratio_AML1_paper];
real_d_final = real_d(selec)./ratio_D6_over_D1;
% real_d_final = real_d(selec);

%% D2
possible_state_label_persite = [[-10;-1;10],[-100;-1;100]];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D2_valid_sublist = valid_sublist;
D2_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% D3
possible_state_label_persite = [[-12;-1;12],[-100;-1;100]];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D3_valid_sublist = valid_sublist;
D3_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% D4
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D4_valid_sublist = valid_sublist;
D4_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% D5
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D5_valid_sublist = valid_sublist;
D5_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% D6 
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D6_valid_sublist = valid_sublist;
D6_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% D7 
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D7_valid_sublist = valid_sublist;
D7_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% F2 
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_12kb = [-12;-1;zeros(length(possible_state_label_14kb)-2,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
F2_valid_sublist = valid_sublist;
F2_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% F3 
possible_state_label_14kb = [-14;-1;3;14;314];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
F3_valid_sublist = valid_sublist;
F3_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% F4 
possible_state_label_14kb = [-14;-1;3;14;314];
possible_state_label_12kb = [-12;-1;zeros(length(possible_state_label_14kb)-2,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
F4_valid_sublist = valid_sublist;
F4_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% E2
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
E2_valid_sublist = valid_sublist;
E2_prob = RNAPHandle(valid_sublist,energy_table,recomputed_concentration_table_Bcell,nuc_destab_location);

%% E3
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
E3_valid_sublist = valid_sublist;
E3_prob = RNAPHandle(valid_sublist,energy_table,recomputed_concentration_table_Bcell,nuc_destab_location);

%% AML1 paper
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
AML1_valid_sublist = valid_sublist;
AML1_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% compute error
sim_d = [D2_prob,D3_prob,D4_prob,D5_prob,D6_prob,D7_prob,F2_prob,F3_prob,F4_prob,E2_prob,E3_prob,AML1_prob];
sim_d_final = sim_d(selec);
% gene_data_err = weighted_msd(sim_d_final,real_d_final);
gene_data_err = 10*weighted_msd(sim_d_final,real_d_final); % "geneDataPenalized"

%% create a dict for better output (I only need D7; save time, not using dict now)
% all_sublist_str.D2 = D2_valid_sublist;
% all_sublist_str.D3 = D3_valid_sublist;
% all_sublist_str.D4 = D4_valid_sublist;
% all_sublist_str.D5 = D5_valid_sublist;
% all_sublist_str.D6 = D6_valid_sublist;
% all_sublist_str.D7 = D7_valid_sublist;
% all_sublist_str.F2 = F2_valid_sublist;
% all_sublist_str.F3 = F3_valid_sublist;
% all_sublist_str.F4 = F4_valid_sublist;
% all_sublist_str.E2 = E2_valid_sublist;
% all_sublist_str.E3 = E3_valid_sublist;

%% enforce D6 and D7
ksmall = D7_prob; 
klarge = D6_prob;

pe = penalty_function(ksmall,klarge,g1,g2);

%% can we find just Hill([pu]), given best fitted pu, to match the best fitted [ce]?
pred_ce_conc = zeros(1,4);
% load('error_prop_0216.mat','remRate_constant_ratio_all_CEBPA')
for i = 1:4
    if i == 3
        remRate = remRate_constant_ratio_all_CEBPA(i,1);
    elseif i == 4
        remRate = 1;
    else
        remRate = remRate_constant_ratio_all_CEBPA(i,2);
    end
    pred_ce_conc(i) = cebpa_RegBy_pu1_Hill(p_cebp_hill,recomputed_pu1_ss_conc(i))/cebpa_RegBy_pu1_Hill(p_cebp_hill,recomputed_pu1_ss_conc(end))/remRate;
end

debug_err = weighted_msd(pred_ce_conc,cebpa_ss_conc);

%% END: make sure production rate with or without mutation is consistent with PU1 gene expression data
%%%%%%%%%%%%%%%%%%%%%%

%% FINAL error!!
% extra penalty: deviation from median also results in error
err1 = weighted_msd(recomputed_pu1_ss_conc(1:3),conc_ratio_all(:,2)');
% extra penalty: deviation from median also results in error
err2 = weighted_msd([recomputed_D7_prob_bcell,recomputed_D7_prob_ProLow,recomputed_D7_prob_ProHigh],prodRate_ratio_all(:,2)');

% if CEPROD == 'Y'
%     final_err = pe + gene_data_err + p1 + p2 + p3 + p4 + p5 + p6 + p21 + p22 + p23 + p24 + p25 + p26 + err1 + err2 + p11 + p12 + p13 + p14 + p15 + p16;
% else
    final_err = pe + gene_data_err + p1 + p2 + p3 + p4 + p5 + p6 + p21 + p22 + p23 + p24 + p25 + p26 + err1 + err2 + debug_err;
% end


