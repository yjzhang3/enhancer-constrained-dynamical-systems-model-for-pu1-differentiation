clear
clc
% model 1 best pars according to bifurcation analysis
A = readmatrix("OneHillPars_fitted_from_noPCfeedback20.txt");
[~,I] = min(A(:,end));
p_cebp_hill = A(I,1:end-1);
p_cebp_hill_Amin = p_cebp_hill(1);
p_cebp_hill_Amax = p_cebp_hill(2);
p_cebp_hill_N = p_cebp_hill(3);
p_cebp_hill_K = p_cebp_hill(4);
MATCONTpars = [p_cebp_hill_Amin,p_cebp_hill_Amax,p_cebp_hill_N,p_cebp_hill_K];
MATCONTpars(1) = -0.0914469;

load('nofeedback_best_results.mat','pu1_ss_conc','cebpa_ss_conc','best_pars')
x0_for_big_model = [pu1_ss_conc;cebpa_ss_conc]';
x0_for_big_model(1,:) = [0.057,0.07];

full_pars = [best_pars(1:12),MATCONTpars];


selec_cebpa = 1:11;
nuc_destab_location_cebpa = [0,1,0,1];

[final_err,pu1_ss_conc,cebpa_ss_conc,pu1_ss_prod,p_cebp_hill,...
    real_d_final,sim_d_final,D7_valid_sublist,energy_table] = ...
ODE_objective_2D_cebpa_PCnoAUTOfeedback_smallestPHdilute_0618_x(selec_cebpa,full_pars,nuc_destab_location_cebpa,x0_for_big_model,'N');