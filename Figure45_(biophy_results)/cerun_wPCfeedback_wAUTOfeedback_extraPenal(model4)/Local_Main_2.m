
%% C and R
clear
clc
close all
selec_cebpa_runx1 = 1:12;
nuc_destab_location_cebpa_runx1 = [0,1,0,1,0,1];

for j = 1:length(nuc_destab_location_cebpa_runx1(:,1))
    % with feedback, with itnersection requirement
    objFuncHandle = @objective_with_cebpa_runx1_general_PU1CEfeedback_0401;
    fit_with_cebpa_runx1_general_PU1CEfeedback(objFuncHandle,nuc_destab_location_cebpa_runx1(j,:),selec_cebpa_runx1); 
end
