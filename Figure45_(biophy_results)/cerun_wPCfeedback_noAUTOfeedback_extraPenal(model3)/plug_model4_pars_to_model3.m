% plug in best-fitted ss for all tf and energies to see whether CEAUTO
% alone can still fit well
clear
clc
load("best_model_ss.mat","energy_table","pu1_ss_conc","cebpa_ss_conc","runx1_ss_conc","best_pars")
selec = 1:12;
nuc_destab_location = [0,1,0,1,0,0];
promoter_only_rescale_factor = best_pars(17);

%% first test using random p_cebp_hill to see if best-fitted energy
% parameters lead to good PU.1 transcription rate prediction and good PU.1 steady state + production rate prediction. 
% YES AGREED. 
p = [0,0,0,0];
[final_err,p_cebp_hill,...
    real_d_final,sim_d_final,D7_valid_sublist] = ...
    ODE_objective_debug_model3_test_only(selec,p,nuc_destab_location,energy_table,promoter_only_rescale_factor,pu1_ss_conc,cebpa_ss_conc,runx1_ss_conc);

%% now use best-fitted ce and pu ss as initial conditions, find p_cebp_hill
% that not only fits to both datasets well like before, but also lead to ce
% ss that agrees with best-fitted cess
clc
lb = [0.0001,1,2,0.001];
ub = [0.1,1,2,2];

fun = @(p) ODE_objective_debug_model3(selec,p,nuc_destab_location,energy_table,promoter_only_rescale_factor,pu1_ss_conc,cebpa_ss_conc,runx1_ss_conc);
fn = 'plug_model4_pars_to_model3_find_p_cebp_hill.txt';
for i=1:10
    [x,fval] = particleswarm(fun,4,lb,ub)
    create_parameter_file(fn,x,fval,'a+');
end