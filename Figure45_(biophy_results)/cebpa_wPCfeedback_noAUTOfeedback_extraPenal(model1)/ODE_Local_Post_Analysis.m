% local post analysis
%% cebpa
delete /outputs_figs
clear
clc
close all 
selec_cebpa = 1:11;
nuc_destab_location_cebpa = [0,1,0,1];

PredFuncHandle = @ODE_prediction_cebpa_PU1CEfeedback;
% ****
objFuncHandle = @ODE_objective_2D_cebpa_PCnoAUTOfeedback_smallestPHdilute_0618;
% objFuncHandle = @ODE_objective_2D_cebpa_PCnoAUTOfeedback_smallestPHdilute_0618_x;

% objFuncHandle = @ODE_objective_2D_cebpa_PU1CEfeedback_smallestPHdilute_0404_v2;
% objFuncHandle = @ODE_objective_2D_with_cebpa_PU1CEfeedback_smallestPHdilute_0401;

% ****
% objFuncHandle = @ODE_objective_with_cebpa_PU1CEfeedback;

% funname = func2str(objFuncHandle);
% orig_name = strcat('ODElocalfit_',funname(15:end),'_PrepC',c'_50N','_maxV1','_UBtau'); 
% fn = strcat(orig_name,'_',num2str(selec_cebpa),'_',num2str(nuc_destab_location_cebpa));
% fn = strrep(fn, ' ', '');
% final_fn = strcat(fn,'.txt');

% ****
final_fn = ...
"ODENOTSfit_PS_fsolve_PCn22_N.txt";   
% 'ODEfit_2D_cebpa_PU1CEfeedback_smallestPHdilute_0404_v2_PattrC_freeN_maxV1_1234567891011_0101.txt';
% 'ODElocalfit_PS_2D_with_cebpa_PU1CEfeedback_smallestPHdilute_0401_PattrC_freeN_maxV1_1234567891011_0101.txt';
% ****this is the final file
CEPROD = final_fn{1}(end-4);

% final_fn = 'ODElocalfit_with_cebpa_PU1CEfeedback_PrepC_10N_maxV1_1234567891011_0101.txt';

for jj = 1:1
    destab_id = num2str(nuc_destab_location_cebpa(jj,:));
    final_id = strrep(destab_id," ","");
    final_selec = strrep(num2str(selec_cebpa)," ","");
    if contains(func2str(objFuncHandle),'x')
        PredFuncHandle(objFuncHandle,selec_cebpa,final_fn,nuc_destab_location_cebpa(jj,:),strcat(CEPROD,'_diffx0'));
    else
    % final_fn = strcat('local_fit_with_cepba_runx1_general_alt1',v_num,final_selec,'_',final_id,'.txt');
    PredFuncHandle(objFuncHandle,selec_cebpa,final_fn,nuc_destab_location_cebpa(jj,:),CEPROD);
    end
    % ODE_predict_cebpa_multiplePars(objFuncHandle,selec_cebpa,final_fn,nuc_destab_location_cebpa(jj,:),CEPROD,6.3);
end


