% % all possible nucleosome mechanisms
% n = 6;
% combinations = dec2bin(0:2^n-1) - '0';
% fl_combinations = flipud(combinations);
% 
% % get slurm task id from environment variable
% task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
% x = fl_combinations(task_id,:);
% 
% % % output path on NOTS
% % NOTS_path = getenv('SHARED_SCRATCH'); % we will store the text files for each nucleosome comboination in $SHARED_SCRATCH
% 
% % run function
% fit_with_cebpa_runx1_general(x)

clear
clc
close all
selec_cebpa_runx1 = 1:12;
nuc_destab_location_cebpa_runx1 = [0,1,0,1,0,0];

%objFuncHandle = @objective_with_cebpa_runx1_general_PU1CEfeedback_0401;
%fit_with_cebpa_runx1_general_PU1CEfeedback(objFuncHandle,nuc_destab_location_cebpa_runx1,selec_cebpa_runx1); 

objFuncHandle = @ODE_objective_2D_with_cebpa_runx1_PU1CEfeedback_smallestPH_0514;
ODE_fit_with_cebpa_runx1_PU1CEfeedback_PS(1,objFuncHandle,nuc_destab_location_cebpa_runx1,selec_cebpa_runx1); 
