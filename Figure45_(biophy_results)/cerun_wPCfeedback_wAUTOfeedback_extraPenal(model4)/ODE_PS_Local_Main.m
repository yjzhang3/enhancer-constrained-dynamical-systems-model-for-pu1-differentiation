%% C and R
clear
clc
close all
selec_cebpa_runx1 = 1:12;
nuc_destab_location_cebpa_runx1 = [0,1,0,1,0,0];

for j = 1:length(nuc_destab_location_cebpa_runx1(:,1))
    % with feedback, with itnersection requirement
    objFuncHandle = @ODE_objective_2D_cerun_PCfeed_CEauto_smallestPH_0618;
    % ODE_fit_with_cebpa_runx1_PU1CEfeedback(objFuncHandle,nuc_destab_location_cebpa_runx1(j,:),selec_cebpa_runx1); 
    ODE_fit_cerun_PCfeed_CEauto_PS(0,objFuncHandle,nuc_destab_location_cebpa_runx1(j,:),selec_cebpa_runx1); 
end
