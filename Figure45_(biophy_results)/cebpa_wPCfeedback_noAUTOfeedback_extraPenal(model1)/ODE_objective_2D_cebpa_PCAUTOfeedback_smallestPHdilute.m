function [final_err,pu1_ss_conc,cebpa_ss_conc,pu1_ss_prod,...
    real_d_final,sim_d_final,D7_valid_sublist,energy_table,...
    Bcell_over_time2,ProLow_over_time2,ProHigh_over_time2,tB2,tPL2,tPH2] = ...
    ODE_objective_2D_cebpa_PCAUTOfeedback_smallestPHdilute(selec,p,nuc_destab_location)
% directly modified from objective_with_cebpa_general_PU1CEfeedback, so see
% information in the other file

%% housekeeping
% load(fullfile(parentdir(pwd, 1), 'Error propagation and Normalization','error_prop_0216.mat')); 
load('error_prop_0216.mat')
load('Cebpa_ct_constraint.mat')

RNAPHandle = @RNAP_occ_prob_cebpa_general;
% D_pu1 = 1/70;
% D_cebpa = log(2)/(100/60); % unit is 1/hour, same as PU1!!

% D7 
site_index = 1;
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D7_valid_sublist = valid_sublist; % need to define early becuase it'll be used in parallel later

%% load pars
G_NUC = p(1);
G_PU1_14 = p(2);
G_PU1_12 = p(3);
G_PU1_10 = p(4);
G_CEBP = p(5);
G_RNAP = p(6);

G_NUC_PU1 = p(7); 
G_NUC_CEBP = p(8);

G_PU1_RNAP = p(9);
G_CEBP_RNAP = p(10);

G_PU1_CEBP = p(11);

promoter_only_rescale_factor = p(12); % this is the ratio between E1 (promoter only in B cell) and D1 (promoter only in macrophage)

% pu1_bcell_conc = p(13); % don't need because i can simulate ODE to get it
% pu1_ProLow_conc = p(14); % I don't need pu1 conc in pro; instead, I can directly simulate ODE to get them, and then constrain it with data
% pu1_ProHigh_conc = p(15);
p_cebp_hill = p(13:16);
p_auto_hill = p(17:20);

energy_table = [G_NUC,G_PU1_14,G_PU1_12,G_PU1_10,G_CEBP,G_RNAP, ... %1-6
    G_NUC_PU1,G_NUC_CEBP,... % 7-8
    G_PU1_RNAP,G_CEBP_RNAP,... % 9-10
    G_PU1_CEBP]; ... % 11
    
% nuc_destab_location =
% [PU1_PRO,PU1_allTFBS,C/EBPa_PRO,C/EBPa_12]


%% construct cell-type specific ODE
% remember D7 is wildtype whose production we're interested in computing
% concentration table order of entry: [nuc],[pu1],[cebp],[RNAP]

ODE_func_Bcell = @(t,y) ODE_func_cellType(t,y,RNAPHandle,p_cebp_hill,p_auto_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(1,2),remRate_constant_ratio_all_CEBPA(1,2));
ODE_func_ProLow = @(t,y) ODE_func_cellType(t,y,RNAPHandle,p_cebp_hill,p_auto_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(2,2),remRate_constant_ratio_all_CEBPA(2,2));
ODE_func_ProHigh = @(t,y) ODE_func_cellType(t,y,RNAPHandle,p_cebp_hill,p_auto_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(3,1),remRate_constant_ratio_all_CEBPA(3,1)); %!!!!!

% event_Bcell = @(t,y) ODE_func_cellType_myEvents(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,tau_Bcell,remRate_constant_ratio_all(1,2),remRate_constant_ratio_all_CEBPA(1,2));
% event_ProLow = @(t,y) ODE_func_cellType_myEvents(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,tau_ProLow,remRate_constant_ratio_all(2,2),remRate_constant_ratio_all_CEBPA(2,2));
% event_ProHigh = @(t,y) ODE_func_cellType_myEvents(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,tau_Bcell,remRate_constant_ratio_all(3,1),remRate_constant_ratio_all_CEBPA(3,1));

%% extract pu1 steady state; precompute ce as a function of pu1 ss; 
% use previously solved pu1ss as initial condition for faster
% performance!!!!
tspan = 200; % original is 200
% tspan = 10;
[tB,Bcell_over_time] = ode23t(ODE_func_Bcell,[0 tspan],[0.01,0.03],odeset('RelTol',1e-10,'AbsTol',1e-12));
[tPL,ProLow_over_time] = ode23t(ODE_func_ProLow,[0 tspan],[0.01,0.04],odeset('RelTol',1e-10,'AbsTol',1e-12));
[tPH,ProHigh_over_time] = ode23t(ODE_func_ProHigh,[0 tspan],[0.43,0.6],odeset('RelTol',1e-10,'AbsTol',1e-12));

%% verify persistence
% eps = 1e-8;
% Tcheck = 400; % original, used for fitting
Tcheck = 1200;
[tB2,Bcell_over_time2] = ode23t(ODE_func_Bcell,[0 Tcheck],Bcell_over_time(end,:)',odeset('RelTol',1e-10,'AbsTol',1e-12));
[tPL2,ProLow_over_time2] = ode23t(ODE_func_ProLow,[0 Tcheck],ProLow_over_time(end,:)',odeset('RelTol',1e-10,'AbsTol',1e-12));
[tPH2,ProHigh_over_time2] = ode23t(ODE_func_ProHigh,[0 Tcheck],ProHigh_over_time(end,:)',odeset('RelTol',1e-10,'AbsTol',1e-12));

% if norm(ProHigh_over_time2(end,:)-ProHigh_over_time(end,:)) < eps
%     disp("PH Stable steady state confirmed")
% else
%     disp("False equilibrium (trajectory still evolving)")
% end
% 
% if norm(ProLow_over_time2(end,:)-ProLow_over_time(end,:)) < eps
%     disp("PL Stable steady state confirmed")
% else
%     disp("False equilibrium (trajectory still evolving)")
% end
% 
% if norm(Bcell_over_time2(end,:)-Bcell_over_time(end,:)) < eps
%     disp("Bcell Stable steady state confirmed")
% else
%     disp("False equilibrium (trajectory still evolving)")
% end

%% construct conc table
pu1_Bcell_ss_conc = Bcell_over_time2(end,1);
pu1_ProLow_ss_conc = ProLow_over_time2(end,1);
pu1_ProHigh_ss_conc = ProHigh_over_time2(end,1);
pu1_ss_conc = [pu1_Bcell_ss_conc,pu1_ProLow_ss_conc,pu1_ProHigh_ss_conc];

cebpa_bcell_ss_conc = Bcell_over_time2(end,2);
cebpa_ProLow_ss_conc = ProLow_over_time2(end,2);
cebpa_ProHigh_ss_conc = ProHigh_over_time2(end,2);
cebpa_ss_conc = [cebpa_bcell_ss_conc,cebpa_ProLow_ss_conc,cebpa_ProHigh_ss_conc];

concentration_table_Mac = [1,1,1,1]; % [nuc],[pu1],[cebp],[RNAP] for *macrophage*
concentration_table_Bcell = [1,pu1_Bcell_ss_conc,cebpa_bcell_ss_conc,1]; % [nuc],[pu1],[cebp],[RNAP]
concentration_table_ProLow = [1,pu1_ProLow_ss_conc,cebpa_ProLow_ss_conc,1];
concentration_table_ProHigh = [1,pu1_ProHigh_ss_conc,cebpa_ProHigh_ss_conc,1]; 

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

%% enforce steady state of CEBPa of each cell type happens at data-constrained range
% Bcell
g1 = 0.05;
g2 = 0.05;
ksmall = Bcell_over_mono_bound(1);
klarge = cebpa_bcell_ss_conc;
p11 = penalty_function(ksmall,klarge,g1,g2);

ksmall = cebpa_bcell_ss_conc;
klarge = Bcell_over_mono_bound(2);
p12 = penalty_function(ksmall,klarge,g1,g2);

% Prolow
g1 = 0.05;
g2 = 0.05;
ksmall = prolow_over_mono_bound(1);
klarge = cebpa_ProLow_ss_conc;
p13 = penalty_function(ksmall,klarge,g1,g2);

ksmall = cebpa_ProLow_ss_conc;
klarge = prolow_over_mono_bound(2);
p14 = penalty_function(ksmall,klarge,g1,g2);

% ProHigh
g1 = 0.05;
g2 = 0.05;
ksmall = prohigh_over_mono_bound(1);
klarge = cebpa_ProHigh_ss_conc;
p15 = penalty_function(ksmall,klarge,g1,g2);

ksmall = cebpa_ProHigh_ss_conc;
klarge = prohigh_over_mono_bound(2);
p16 = penalty_function(ksmall,klarge,g1,g2);

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

% ratio_AML1_paper = ratio_D4_over_D1*32/880; %12 % without AML1 binding sites, decrease in PU1 expression
% ratio_CEBP_paper_C1mut = 77/100; % do you really need all CEBPa binding sites? The effect of taking out all 3 cebap sites and taking out the whole -14kb region is the same!
% ratio_CEBP_paper_C2mut = 59/100;
% ratio_CEBP_paper_C5mut = 20/100;

real_d = [ratio_D2_over_D1,ratio_D3_over_D1,ratio_D4_over_D1,ratio_D5_over_D1,ratio_D6_over_D1,ratio_D7_over_D1,ratio_F2_over_D1,ratio_F3_over_D1,ratio_F4_over_D1,ratio_E2_over_D1,ratio_E3_over_D1];
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
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D4_valid_sublist = valid_sublist;
D4_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% D5
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D5_valid_sublist = valid_sublist;
D5_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% D6 
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D6_valid_sublist = valid_sublist;
D6_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% D7 
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D7_valid_sublist = valid_sublist;
D7_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% F2 
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_12kb = [-12;-1;zeros(length(possible_state_label_14kb)-2,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
F2_valid_sublist = valid_sublist;
F2_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% F3 
possible_state_label_14kb = [-14;-1;14];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
F3_valid_sublist = valid_sublist;
F3_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% F4 
possible_state_label_14kb = [-14;-1;14];
possible_state_label_12kb = [-12;-1;zeros(length(possible_state_label_14kb)-2,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
F4_valid_sublist = valid_sublist;
F4_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Mac,nuc_destab_location);

%% E2 (B cell)
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
E2_valid_sublist = valid_sublist;
E2_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Bcell,nuc_destab_location);

%% E3 (B cell)
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
E3_valid_sublist = valid_sublist;
E3_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Bcell,nuc_destab_location);

%% compute error
sim_d = [D2_prob,D3_prob,D4_prob,D5_prob,D6_prob,D7_prob,F2_prob,F3_prob,F4_prob,E2_prob,E3_prob];
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
final_err = gene_data_err + p1 + p2 + p3 + p4 + p5 + p6 + p11 + p12 + p13 + p14 + p15 + p16 + p21 + p22 + p23 + p24 + p25 + p26;

% %%
% %%% nested function
% function status = outFcn(t, y, flag)
%     status = 0;   % never stop solver manually
%     if strcmp(flag,'done')
%         y_final = y(end);
%         y2_at_stop = y_final^2;
%     end
% end
