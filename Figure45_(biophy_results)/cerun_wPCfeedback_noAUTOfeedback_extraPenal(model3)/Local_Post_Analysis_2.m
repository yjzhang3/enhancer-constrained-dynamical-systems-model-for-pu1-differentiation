% local post analysis
%% cebpa and runx1
clear
close all
clc
selec_cebpa_runx1 = 1:12;
nuc_destab_location_cebpa_runx1 = [0,1,0,1,0,0];

PredFuncHandle = @predict_cebpa_runx1_PU1CEfeedback;
TestHillFuncHandle = @test_Hill_property_cebpa_runx1_general_PU1CEfeedback;
objFuncHandle = @objective_with_cebpa_runx1_general_PU1CEfeedback_0401;

funname = func2str(objFuncHandle);
orig_name = strcat('local_fit_',funname(11:end),'_PattrC','_2N','_maxV1');  
fn = strcat(orig_name,'_',num2str(selec_cebpa_runx1),'_',num2str(nuc_destab_location_cebpa_runx1));
fn = strrep(fn, ' ', '');
final_fn = strcat(fn,'.txt');

for jj = 1:1
    destab_id = num2str(nuc_destab_location_cebpa_runx1(jj,:));
    final_id = strrep(destab_id," ","");
    final_selec = strrep(num2str(selec_cebpa_runx1)," ","");
    % final_fn = strcat('local_fit_with_cepba_runx1_general_alt1',v_num,final_selec,'_',final_id,'.txt');
    final_fn1 = strcat('outputs_mats/',fn,'.mat');
    PredFuncHandle(objFuncHandle,selec_cebpa_runx1,final_fn,nuc_destab_location_cebpa_runx1(jj,:));
    TestHillFuncHandle(objFuncHandle,final_fn1,nuc_destab_location_cebpa_runx1(jj,:));
end