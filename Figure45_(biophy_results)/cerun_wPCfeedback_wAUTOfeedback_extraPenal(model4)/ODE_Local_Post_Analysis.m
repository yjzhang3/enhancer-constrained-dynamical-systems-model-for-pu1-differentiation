% local post analysis
%% cebpa and runx1
clear
close all
clc
selec_cebpa_runx1 = 1:12;
nuc_destab_location_cebpa_runx1 = [0,1,0,1,0,0];

PredFuncHandle = @ODE_prediction_cebpa_runx1_PU1CEfeedback;
% ****close all

objFuncHandle = @ODE_objective_2D_cerun_PCfeed_CEauto_smallestPH_0618;
% ****

% funname = func2str(objFuncHandle);
% orig_name = strcat('ODElocalfit_PS_',funname(15:end),'_smallestPLdilute_','PrepC','_2N','_maxV1'); 
% fn = strcat(orig_name,'_',num2str(selec_cebpa_runx1),'_',num2str(nuc_destab_location_cebpa_runx1));
% fn = strrep(fn, ' ', '');
% final_fn = strcat(fn,'.txt');
% ****
final_fn = "ODENOTSfit_PS_fsolve_extraP_PCn22_AUTOn44_RPRC00_N.txt";
CEPROD = final_fn{1}(end-4);
% final_fn = 'ODENOTSfit_PS_10.txt';
% ****

for jj = 1:1
    destab_id = num2str(nuc_destab_location_cebpa_runx1(jj,:));
    final_id = strrep(destab_id," ","");
    final_selec = strrep(num2str(selec_cebpa_runx1)," ","");
    PredFuncHandle(objFuncHandle,selec_cebpa_runx1,final_fn,nuc_destab_location_cebpa_runx1(jj,:),CEPROD);
    % ODE_predict_cerun_multiplePars(objFuncHandle,selec_cebpa_runx1,final_fn,nuc_destab_location_cebpa_runx1(jj,:),CEPROD,0.32)
end