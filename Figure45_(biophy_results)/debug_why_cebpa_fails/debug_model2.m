% to debug why the presence of RUNX1 must be necessary, let's fix energy
% pars in pu1 production (using model 2 best pars, because 
% gene expressio is fitted well in model2), as well as using pu1 levels that are
% consistent with experiment (adjust the pro high level, because model 2
% can't predict correc prohigh pu1 level), and fit for cebpa LEVELS
clear
clc
load("model2_best_results.mat",'best_pars','pu1_ss_conc','energy_table',"cebpa_ss_conc")
promoter_only_rescale_factor = best_pars(12);
selec = 1:11;
nuc_destab_location = [0,1,0,1];

%% first see if energy pars from model 2 lead to good gene expression
% prediction (to do this, you also need the originally fitted cebpa
% concentration!). Confirmed it was correct
ODE_objective_debug_model2_test_only(selec,[],nuc_destab_location,pu1_ss_conc,cebpa_ss_conc,energy_table,promoter_only_rescale_factor)

%% now slightly adjust pu1_ss_conc, and solve for cebpa_ss, by fitting p_auto_hill and p_cepba_hill
