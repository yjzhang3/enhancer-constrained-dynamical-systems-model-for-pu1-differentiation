function [final_er,sim_d_final,real_d_final,energy_table,...
    concentration_table,concentration_table_Bcell,...
    concentration_table_ProLow,concentration_table_ProHigh,all_sublist_str,...
    norm_D7_prob_bcell,norm_D7_prob_ProLow,norm_D7_prob_ProHigh,p_cebp_hill,pu1_prod_slope_celltype] = ...
    objective_cebpa_general_PU1CEfeedback_smallestPHdilute_0401(selec,p,changePU1conc,pu1Macconc,pu1Bconc,nuc_destab_location)

%% now add back in the itnersection requirment
%%
load(fullfile(parentdir(pwd, 1), 'Error propagation and Normalization', 'error_prop_0216.mat'));
RNAPHandle = @RNAP_occ_prob_cebpa_general;
% in this model, 
% C/EBPa, and PU1 are TF
% C/EBPa are dimers!!

% we will have all possibitilies of transcriptional mechanism I can think
% of, and then set these parameters to zero in the fitting function. 
% * PU1 destablizes nucleosome at all PU.1 binding sites
% * PU1 destablizes nucleosome at promoter
% * PU1 attracts RNAP, regardless of which PU1 site
% * C/EBPa destabilizes nucleosome bound at -12kb only
% * C/EBPa destabilizes nucleosome bound at the promoter
% * C/EBPa recruits RNAP


% C/EBPa repels
% PU1
%% unknown pars
% -100 = nucleosome bound at promoter
% -14, -12, -10 = nucleosome bound at TFBS -14, -12, -10
% -1 = anysite unbound
% 12 = PU.1 bound at -12
% 10 = PU.1 bound at -10
% 100 = RNAP bound

% now enumrate all possibilities for the -14 (in embryotic stem cells,
% RUNX1 expresses first, then PU1, then C/EBPa)
% 2 = CEBPA bound alone

% 14 = PU.1 bound alone at -14kb
% 214 = PU.1 and C/EBPa bound together


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

promoter_only_rescale_factor = p(12);

pu1_bcell_conc = p(13);
pu1_ProLow_conc = p(14);
pu1_ProHigh_conc = p(15);
p_cebp_hill = p(16:19);

energy_table = [G_NUC,G_PU1_14,G_PU1_12,G_PU1_10,G_CEBP,G_RNAP, ... %1-6
    G_NUC_PU1,G_NUC_CEBP,... % 7-8
    G_PU1_RNAP,G_CEBP_RNAP,... % 9-10
    G_PU1_CEBP]; ... % 11
    
% nuc_destab_location =
% [PU1_PRO,PU1_allTFBS,C/EBPa_PRO,C/EBPa_12]

% precompute cebpa concentration needed for b cell, ProHigh and ProLow
D_cebpa = log(2)/(100/60); % unit is 1/hour, same as PU1!!
tau_Bcell = 25.5;
tau_Mac = 40;
tau_ProLow = 28.8;
tau_ProHigh = 11.6;

cebpa_bcell_conc = (cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_bcell_conc)/cebpa_RegBy_pu1_Hill(p_cebp_hill,1))/...
    remRate_constant_ratio_all_CEBPA(1,2);
cebpa_ProLow_conc = (cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_ProLow_conc)/cebpa_RegBy_pu1_Hill(p_cebp_hill,1))/...
    remRate_constant_ratio_all_CEBPA(2,2);
cebpa_ProHigh_conc = (cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_ProHigh_conc)/cebpa_RegBy_pu1_Hill(p_cebp_hill,1))/...
    remRate_constant_ratio_all_CEBPA(3,1);

if changePU1conc == 0
    concentration_table = [1,1,1,1]; % [nuc],[pu1],[cebp],[runx],[RNAP] for *macrophage*
    concentration_table_Bcell = [1,pu1_bcell_conc,cebpa_bcell_conc,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
elseif changePU1conc == 1 % when I wanto change PU1 conc in mac to look at TR dependence on PU1 conc 
    concentration_table = [1,pu1Macconc,1,1]; % [nuc],[pu1],[cebp],[runx],[RNAP] for *macrophage*
    concentration_table_Bcell = [1,pu1_bcell_conc,cebpa_bcell_conc,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
elseif changePU1conc == 2 % when I wanto change PU1 conc in b cell to look at TR dependence on PU1 conc 
    concentration_table = [1,1,1,1]; % for *macrophage*
    concentration_table_Bcell = [1,pu1Bconc,cebpa_bcell_conc,1]; % [nuc],[pu1],[cebp],[RNAP]
end
site_index = 1;

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
D2_prob = RNAPHandle(valid_sublist,energy_table,concentration_table,nuc_destab_location);

%% D3
possible_state_label_persite = [[-12;-1;12],[-100;-1;100]];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D3_valid_sublist = valid_sublist;
D3_prob = RNAPHandle(valid_sublist,energy_table,concentration_table,nuc_destab_location);

%% D4
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D4_valid_sublist = valid_sublist;
D4_prob = RNAPHandle(valid_sublist,energy_table,concentration_table,nuc_destab_location);

%% D5
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D5_valid_sublist = valid_sublist;
D5_prob = RNAPHandle(valid_sublist,energy_table,concentration_table,nuc_destab_location);

%% D6 
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D6_valid_sublist = valid_sublist;
D6_prob = RNAPHandle(valid_sublist,energy_table,concentration_table,nuc_destab_location);

%% D7 
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D7_valid_sublist = valid_sublist;
D7_prob = RNAPHandle(valid_sublist,energy_table,concentration_table,nuc_destab_location);

%% F2 
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_12kb = [-12;-1;zeros(length(possible_state_label_14kb)-2,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
F2_valid_sublist = valid_sublist;
F2_prob = RNAPHandle(valid_sublist,energy_table,concentration_table,nuc_destab_location);

%% F3 
possible_state_label_14kb = [-14;-1;14];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
F3_valid_sublist = valid_sublist;
F3_prob = RNAPHandle(valid_sublist,energy_table,concentration_table,nuc_destab_location);

%% F4 
possible_state_label_14kb = [-14;-1;14];
possible_state_label_12kb = [-12;-1;zeros(length(possible_state_label_14kb)-2,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
F4_valid_sublist = valid_sublist;
F4_prob = RNAPHandle(valid_sublist,energy_table,concentration_table,nuc_destab_location);

%% E2
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
E2_valid_sublist = valid_sublist;
E2_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Bcell,nuc_destab_location);

%% E3
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
E3_valid_sublist = valid_sublist;
E3_prob = RNAPHandle(valid_sublist,energy_table,concentration_table_Bcell,nuc_destab_location);

%% AML1 paper
% possible_state_label_14kb = [-14;-1;2;14;214];
% possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
% possible_state_label_persite = [possible_state_label_14kb,possible_state_label_PRO];
% valid_sublist = sublist_finder(possible_state_label_persite,site_index);
% AML1_valid_sublist = valid_sublist;
% AML1_prob = RNAPHandle(valid_sublist,energy_table,concentration_table,nuc_destab_location);

%% all simulated data 
% sim_d = exp(coeff)*[D2_prob,D3_prob,D4_prob,D5_prob,D6_prob,D7_prob,F2_prob,F3_prob,F4_prob,E2_prob,E3_prob];
sim_d = [D2_prob,D3_prob,D4_prob,D5_prob,D6_prob,D7_prob,F2_prob,F3_prob,F4_prob,E2_prob,E3_prob];
sim_d_final = sim_d(selec);

%% error
% er = weighted_msd(sim_d_final,real_d_final);
er = 10*weighted_msd(sim_d_final,real_d_final); % "geneDataPenalized"

%% create a dict for better output
all_sublist_str.D2 = D2_valid_sublist;
all_sublist_str.D3 = D3_valid_sublist;
all_sublist_str.D4 = D4_valid_sublist;
all_sublist_str.D5 = D5_valid_sublist;
all_sublist_str.D6 = D6_valid_sublist;
all_sublist_str.D7 = D7_valid_sublist;
all_sublist_str.F2 = F2_valid_sublist;
all_sublist_str.F3 = F3_valid_sublist;
all_sublist_str.F4 = F4_valid_sublist;
all_sublist_str.E2 = E2_valid_sublist;
all_sublist_str.E3 = E3_valid_sublist;
% all_sublist_str.AML1 = AML1_valid_sublist;

%% penalty 1
ksmall = D7_prob; 
klarge = D6_prob;
g1 = 0.05;
g2 = 0.05;

penalty_weight1 = penalty_function(ksmall,klarge,g1,g2);

%% penalty 2
g1 = 0.05;
g2 = 0.05;
% ksmall = abs(G_RUNX_PU1); % less negative (but make sure both are G negative)
% klarge = abs(G_RUNX_CEBP); % more negative
% penalty_weight2 = penalty_function(ksmall,klarge,g1,g2);
penalty_weight2 = 0;

%% penalty 3 4 5 for B cell
D = 1/70;
tau_Bcell = 25.5;
tau_Mac = 40;

norm_remv_rate_at_bcell_pu1 = remRate_constant_ratio_all(1,2)*pu1_bcell_conc; % normalized!

D7_prob_bcell = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_Bcell,nuc_destab_location);

norm_D7_prob_bcell = D7_prob_bcell/D7_prob;

% *normalized* intersect at pu1 bcell
% penalty_weight3 = (norm_D7_prob_bcell-norm_growth_rate_at_bcell_pu1)^2;
penalty_weight3 = 10*(norm_D7_prob_bcell-norm_remv_rate_at_bcell_pu1)^2/(norm_D7_prob_bcell)^2;

% lb < bcell/mac pu1 normalized production rate ratio < ub
g1 = 0.05;
g2 = 0.05;
ksmall = prodRate_ratio_all(1,1);
klarge = norm_D7_prob_bcell;
penalty_weight4 = penalty_function(ksmall,klarge,g1,g2);

ksmall = norm_D7_prob_bcell;
klarge = prodRate_ratio_all(1,3);
penalty_weight5 = penalty_function(ksmall,klarge,g1,g2);

%% penalty 6-11 for pro low and pro high
% non-PU1 TF concentration in progenitor (normalized with respect to 1, the
% concentration in macrophage)
tau_ProLow = 28.8;
% tau_ProLow_Q1 = 19.9;
% tau_ProLow_Q3 = 47.8;
tau_ProHigh = 11.6;
% tau_ProHigh_Q1 = 9.47;
% tau_ProHigh_Q3 = 15.8;

% normalized growth rate of Pro Low
norm_remv_rate_at_ProLow_pu1 = remRate_constant_ratio_all(2,2)*pu1_ProLow_conc;


% normalized growth rate of Pro High 
norm_remv_rate_at_ProHigh_pu1 = remRate_constant_ratio_all(3,1)*pu1_ProHigh_conc; % normalized!


% normalized production rate of Pro Low 
concentration_table_ProLow = [1,pu1_ProLow_conc,cebpa_ProLow_conc,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
D7_prob_ProLow = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_ProLow,nuc_destab_location);
norm_D7_prob_ProLow = D7_prob_ProLow/D7_prob;

% normalized production rate of Pro Low 
concentration_table_ProHigh = [1,pu1_ProHigh_conc,cebpa_ProHigh_conc,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
D7_prob_ProHigh = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_ProHigh,nuc_destab_location);
norm_D7_prob_ProHigh = D7_prob_ProHigh/D7_prob;

% directly let growth rate = production rate
penalty_weight6 = 10*(norm_D7_prob_ProLow-norm_remv_rate_at_ProLow_pu1)^2/(norm_D7_prob_ProLow)^2;
penalty_weight7 = 10*(norm_D7_prob_ProHigh-norm_remv_rate_at_ProHigh_pu1)^2/(norm_D7_prob_ProHigh)^2;

% lb < ProLow/mac pu1 normalized production rate ratio < ub
g1 = 0.05;
g2 = 0.05;
ksmall = prodRate_ratio_all(2,1);
klarge = norm_D7_prob_ProLow;
penalty_weight8 = penalty_function(ksmall,klarge,g1,g2);

ksmall = norm_D7_prob_ProLow;
klarge = prodRate_ratio_all(2,3);
penalty_weight9 = penalty_function(ksmall,klarge,g1,g2);

% lb < ProHigh/mac pu1 normalized production rate ratio < ub
g1 = 0.05;
g2 = 0.05;
ksmall = prodRate_ratio_all(3,1);
klarge = norm_D7_prob_ProHigh;
penalty_weight10 = penalty_function(ksmall,klarge,g1,g2);

ksmall = norm_D7_prob_ProHigh;
klarge = prodRate_ratio_all(3,3);
penalty_weight11 = penalty_function(ksmall,klarge,g1,g2);

%% penatly 12-14 (ensure intersection is stable)
delta = 0.0001;

% growth rate slope of Prolow and ProHigh
norm_remv_rate_at_bcell_pu1_slope = remRate_constant_ratio_all(1,2); 
norm_remv_rate_at_ProLow_pu1_slope = remRate_constant_ratio_all(2,2); 
norm_remv_rate_at_ProHigh_pu1_slope = remRate_constant_ratio_all(3,1); % normalized!

% estimated slope of production rate at Bcell concentration
concentration_table_bcell_left = [1,pu1_bcell_conc-delta,cebpa_bcell_conc,1]; % [nuc],[pu1],[cebp],[RNAP]
D7_prob_bcell_left = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_bcell_left,nuc_destab_location);
norm_D7_prob_bcell_left = D7_prob_bcell_left/D7_prob;

concentration_table_bcell_right = [1,pu1_bcell_conc+delta,cebpa_bcell_conc,1]; % [nuc],[pu1],[cebp],[RNAP]
D7_prob_bcell_right = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_bcell_right,nuc_destab_location);
norm_D7_prob_bcell_right = D7_prob_bcell_right/D7_prob;

norm_D7_prob_bcell_slope = (norm_D7_prob_bcell_right-norm_D7_prob_bcell_left)/(2*delta);

% estimated slope of production rate at ProLow concentration
concentration_table_ProLow_left = [1,pu1_ProLow_conc-delta,cebpa_ProLow_conc,1]; % [nuc],[pu1],[cebp],[RNAP]
D7_prob_ProLow_left = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_ProLow_left,nuc_destab_location);
norm_D7_prob_ProLow_left = D7_prob_ProLow_left/D7_prob;

concentration_table_ProLow_right = [1,pu1_ProLow_conc+delta,cebpa_ProLow_conc,1]; % [nuc],[pu1],[cebp],[RNAP]
D7_prob_ProLow_right = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_ProLow_right,nuc_destab_location);
norm_D7_prob_ProLow_right = D7_prob_ProLow_right/D7_prob;

norm_D7_prob_ProLow_slope = (norm_D7_prob_ProLow_right-norm_D7_prob_ProLow_left)/(2*delta);

% estimated slope of production rate at ProHigh concentration
concentration_table_ProHigh_left = [1,pu1_ProHigh_conc-delta,cebpa_ProHigh_conc,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
D7_prob_ProHigh_left = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_ProHigh_left,nuc_destab_location);
norm_D7_prob_ProHigh_left = D7_prob_ProHigh_left/D7_prob;

concentration_table_ProHigh_right = [1,pu1_ProHigh_conc+delta,cebpa_ProHigh_conc,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
D7_prob_ProHigh_right = RNAPHandle(D7_valid_sublist,energy_table,concentration_table_ProHigh_right,nuc_destab_location);
norm_D7_prob_ProHigh_right = D7_prob_ProHigh_right/D7_prob;

norm_D7_prob_ProHigh_slope = (norm_D7_prob_ProHigh_right-norm_D7_prob_ProHigh_left)/(2*delta);

% bcell production sleop < growth rate slope at intersection
g1 = 0.05;
g2 = 0.05;
ksmall = norm_D7_prob_bcell_slope;
klarge = norm_remv_rate_at_bcell_pu1_slope;
penalty_weight12 = penalty_function(ksmall,klarge,g1,g2);

% ProLow production rate slope < growth rate slope at intersection
g1 = 0.05;
g2 = 0.05;
ksmall = norm_D7_prob_ProLow_slope;
klarge = norm_remv_rate_at_ProLow_pu1_slope;
penalty_weight13 = penalty_function(ksmall,klarge,g1,g2);

% ProHigh production rate slope < growth rate slope at intersection
ksmall = norm_D7_prob_ProHigh_slope;
klarge = norm_remv_rate_at_ProHigh_pu1_slope;
penalty_weight14 = penalty_function(ksmall,klarge,g1,g2);

pu1_prod_slope_celltype = [norm_D7_prob_bcell_slope,norm_D7_prob_ProLow_slope,norm_D7_prob_ProHigh_slope];

%% final error

final_er = er+penalty_weight1+penalty_weight2+penalty_weight3+penalty_weight4+penalty_weight5+penalty_weight6+penalty_weight7+penalty_weight8+penalty_weight9+penalty_weight10+penalty_weight11+penalty_weight12+penalty_weight13+penalty_weight14;
end 
