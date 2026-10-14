%% C only
clear
clc
close all
selec_cebpa = 1:11;
% nuc_destab_location = [[0,1,0,1,0,0];[0,1,0,1,0,1]]; % when there is
% RUNX1 in the model
nuc_destab_location_cebpa = [0,1,0,1]; % there is no RUNX1 in the model

for j = 1:length(nuc_destab_location_cebpa(:,1))
    % no feedback, no itnersection requirement
    % objFuncHandle = @objective_with_cebpa_general_noPCfeedback_geneonly;
    % fit_with_cebpa_general_noPCfeedback_geneonly(objFuncHandle,nuc_destab_location_cebpa(j,:),selec_cebpa); 

    % with feedback, no itnersection requirement
    % objFuncHandle = @objective_with_cebpa_general_PU1CEfeedback_geneonly;
    % fit_with_cebpa_general_PU1CEfeedback_geneonly(objFuncHandle,nuc_destab_location_cebpa(j,:),selec_cebpa); 

    % no feedback, intersection requirement
    % objFuncHandle = @objective_with_cebpa_general_noPCfeedback;
    % fit_with_cebpa_general_noPCfeedback(objFuncHandle,nuc_destab_location_cebpa(j,:),selec_cebpa); 

    % with feedback, with itnersection requirement
    % objFuncHandle = @objective_with_cebpa_general_PU1CEfeedback;
    % fit_with_cebpa_general_PU1CEfeedback(objFuncHandle,nuc_destab_location_cebpa(j,:),selec_cebpa); 

    % with feedback, with itnersection requirement, use ODE simulation
%     objFuncHandle = @ODE_objective_with_cebpa_PU1CEfeedback;
    objFuncHandle = @ODE_objective_2D_cebpa_PU1CEfeedback_smallestPHdilute_0505;
    ODE_fit_PS_with_cebpa_general_PU1CEfeedback(objFuncHandle,nuc_destab_location_cebpa(j,:),selec_cebpa);
end

