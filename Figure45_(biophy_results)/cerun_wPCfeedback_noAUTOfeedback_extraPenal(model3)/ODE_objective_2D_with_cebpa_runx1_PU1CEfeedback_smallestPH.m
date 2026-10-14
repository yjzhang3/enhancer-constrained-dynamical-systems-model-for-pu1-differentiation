function [final_err,pu1_ss_conc,cebpa_ss_conc,runx1_ss_conc,pu1_ss_prod,...
    real_d_final,sim_d_final,D7_valid_sublist,energy_table,...
    Bcell_conc_over_time,ProLow_conc_over_time,ProHigh_conc_over_time,tB,tPL,tPH] = ...
    ODE_objective_2D_with_cebpa_runx1_PU1CEfeedback_smallestPH(selec,p,nuc_destab_location,conc_ratio_all,prodRate_ratio_all)
% directly modified from objective_with_cebpa_general_PU1CEfeedback, so see
% information in the other file

%% housekeeping
load(fullfile(parentdir(pwd, 1), 'Error propagation and Normalization','error_prop_0216.mat')); 
RNAPHandle = @RNAP_occ_prob_cebpa_runx1_general;
D_pu1 = 1/70;
D_cebpa = log(2)/(100/60); % unit is 1/hour, same as PU1!!
tau_Bcell = 25.5;
tau_Mac = 40;
tau_ProLow = 28.8;
tau_ProHigh = 11.6;
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
G_NUC = p(1);
G_PU1_14 = p(2);
G_PU1_12 = p(3);
G_PU1_10 = p(4);
G_CEBP = p(5);
G_RUNX = p(6);
G_RNAP = p(7);
G_NUC_PU1 = p(8); 
G_NUC_CEBP = p(9);
G_NUC_RUNX1 = p(10);
G_PU1_RNAP = p(11);
G_CEBP_RNAP = p(12);
G_RUNX1_RNAP = p(13);
G_RUNX_PU1 = p(14);
G_RUNX_CEBP = p(15);
G_PU1_CEBP = p(16);

promoter_only_rescale_factor = p(17);
runx1_bcell_conc = p(18);
runx1_ProLow_conc = p(19);
runx1_ProHigh_conc = p(20);

p_cebp_hill = p(21:24);

energy_table = [G_NUC,G_PU1_14,G_PU1_12,G_PU1_10,G_CEBP,G_RUNX,G_RNAP, ... %1-7
    G_NUC_PU1,G_NUC_CEBP,G_NUC_RUNX1,... % 8-10
    G_PU1_RNAP,G_CEBP_RNAP,G_RUNX1_RNAP,... % 11-13
    G_RUNX_PU1,G_RUNX_CEBP,G_PU1_CEBP]; ... % 14-16

% nuc_destab_location =
% [PU1_PRO,PU1_allTFBS,C/EBPa_PRO,C/EBPa_12,RUNX1_PRO,RUNX1_12]

%% construct cell-type specific ODE
% remember D7 is wildtype whose production we're interested in computing
% concentration table order of entry: [nuc],[pu1],[cebp],[runx],[RNAP]

ODE_func_Bcell = @(t,y) ODE_func_cellType(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,tau_Bcell,remRate_constant_ratio_all(1,2),remRate_constant_ratio_all_CEBPA(1,2),runx1_bcell_conc);
ODE_func_ProLow = @(t,y) ODE_func_cellType(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,tau_ProLow,remRate_constant_ratio_all(2,2),remRate_constant_ratio_all_CEBPA(2,2),runx1_ProLow_conc);
ODE_func_ProHigh = @(t,y) ODE_func_cellType(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,tau_ProHigh,remRate_constant_ratio_all(3,1),remRate_constant_ratio_all_CEBPA(3,1),runx1_ProHigh_conc);

%% extract pu1 steady state; precompute ce as a function of pu1 ss; 
% use previously solved pu1ss as initial condition for faster
% performance!!!!
tspan = 400;
[tB,Bcell_conc_over_time] = ode45(ODE_func_Bcell,[0 tspan],[0.017,0.001],odeset('RelTol',1e-5));
[tPL,ProLow_conc_over_time] = ode45(ODE_func_ProLow,[0 tspan],[0.18,0.036],odeset('RelTol',1e-5));
[tPH,ProHigh_conc_over_time] = ode45(ODE_func_ProHigh,[0 tspan],[0.373,0.214],odeset('RelTol',1e-5));

pu1_Bcell_ss_conc = Bcell_conc_over_time(end,1);
pu1_ProLow_ss_conc = ProLow_conc_over_time(end,1);
pu1_ProHigh_ss_conc = ProHigh_conc_over_time(end,1);
pu1_ss_conc = [pu1_Bcell_ss_conc,pu1_ProLow_ss_conc,pu1_ProHigh_ss_conc];

cebpa_bcell_ss_conc = Bcell_conc_over_time(end,2);
cebpa_ProLow_ss_conc = ProLow_conc_over_time(end,2);
cebpa_ProHigh_ss_conc = ProHigh_conc_over_time(end,2);
cebpa_ss_conc = [cebpa_bcell_ss_conc,cebpa_ProLow_ss_conc,cebpa_ProHigh_ss_conc];

concentration_table_Mac = [1,1,1,1,1]; % [nuc],[pu1],[cebp],[runx1],[RNAP] for *macrophage*
concentration_table_Bcell = [1,pu1_Bcell_ss_conc,cebpa_bcell_ss_conc,runx1_bcell_conc,1]; 
concentration_table_ProLow = [1,pu1_ProLow_ss_conc,cebpa_ProLow_ss_conc,runx1_ProLow_conc,1];
concentration_table_ProHigh = [1,pu1_ProHigh_ss_conc,cebpa_ProHigh_ss_conc,runx1_ProHigh_conc,1]; 

%% enforce PU1 to be smallest to largest from b cell to pro low, pro high (come on, this is alaready enforced below when compared to real data!)
% g1 = 1;
% g2 = 1;
% ksmall = pu1_Bcell_ss_conc;
% klarge = pu1_ProLow_ss_conc;
% pa = penalty_function(ksmall,klarge,g1,g2);
% 
% ksmall = pu1_ProLow_ss_conc;
% klarge = pu1_ProHigh_ss_conc;
% pb = penalty_function(ksmall,klarge,g1,g2);

%% enforce steady state of PU1 of each cell type happens at data-constrained range
% Bcell
g1 = 0.05;
g2 = 0.05;
ksmall = conc_ratio_all(1,1);
klarge = pu1_Bcell_ss_conc;
p1 = penalty_function(ksmall,klarge,g1,g2);

ksmall = pu1_Bcell_ss_conc;
klarge = conc_ratio_all(1,3);
p2 = penalty_function(ksmall,klarge,g1,g2);

% Prolow
g1 = 0.05;
g2 = 0.05;
ksmall = conc_ratio_all(2,1);
klarge = pu1_ProLow_ss_conc;
p3 = penalty_function(ksmall,klarge,g1,g2);

ksmall = pu1_ProLow_ss_conc;
klarge = conc_ratio_all(2,3);
p4 = penalty_function(ksmall,klarge,g1,g2);

% ProHigh
g1 = 0.05;
g2 = 0.05;
ksmall = conc_ratio_all(3,1);
klarge = pu1_ProHigh_ss_conc;
p5 = penalty_function(ksmall,klarge,g1,g2);

ksmall = pu1_ProHigh_ss_conc;
klarge = conc_ratio_all(3,3);
p6 = penalty_function(ksmall,klarge,g1,g2);

%% enforce production lies in data range
% Bcell
D7_prob_bcell = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_Bcell,nuc_destab_location)/...
    RNAPHandle(D7_valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);
g1 = 0.05;
g2 = 0.05;
ksmall = prodRate_ratio_all(1,1);
klarge = D7_prob_bcell;
p21 = penalty_function(ksmall,klarge,g1,g2);

ksmall = D7_prob_bcell;
klarge = prodRate_ratio_all(1,3);
p22 = penalty_function(ksmall,klarge,g1,g2);

% Prolow
D7_prob_ProLow = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_ProLow,nuc_destab_location)/...
    RNAPHandle(D7_valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);
g1 = 0.05;
g2 = 0.05;
ksmall = prodRate_ratio_all(2,1);
klarge = D7_prob_ProLow;
p23 = penalty_function(ksmall,klarge,g1,g2);

ksmall = D7_prob_ProLow;
klarge = prodRate_ratio_all(2,3);
p24 = penalty_function(ksmall,klarge,g1,g2);

% ProHigh
D7_prob_ProHigh = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_ProHigh,nuc_destab_location)/...
    RNAPHandle(D7_valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);
g1 = 0.05;
g2 = 0.05;
ksmall = prodRate_ratio_all(3,1);
klarge = D7_prob_ProHigh;
p25 = penalty_function(ksmall,klarge,g1,g2);

ksmall = D7_prob_ProHigh;
klarge = prodRate_ratio_all(3,3);
p26 = penalty_function(ksmall,klarge,g1,g2);

pu1_ss_prod = [D7_prob_bcell,D7_prob_ProLow,D7_prob_ProHigh];
runx1_ss_conc = p(18:20);
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
E2_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Bcell,nuc_destab_location);

%% E3
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
E3_valid_sublist = valid_sublist;
E3_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Bcell,nuc_destab_location);

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

%% END: make sure production rate with or without mutation is consistent with PU1 gene expression data
%%%%%%%%%%%%%%%%%%%%%%

%% FINAL error!!
final_err = gene_data_err + p1 + p2 + p3 + p4 + p5 + p6 + p21 + p22 + p23 + p24 + p25 + p26;



