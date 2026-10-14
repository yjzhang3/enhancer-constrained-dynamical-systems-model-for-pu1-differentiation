function ODE_PS_NOTS_Main(moreTime,PCn,CEPROD)
selec_cebpa_runx1 = 1:12;
nuc_destab_location_cebpa_runx1 = [0,1,0,1,0,0];

%objFuncHandle = @objective_with_cebpa_runx1_general_PU1CEfeedback_0401;
%fit_with_cebpa_runx1_general_PU1CEfeedback(objFuncHandle,nuc_destab_location_cebpa_runx1,selec_cebpa_runx1); 
% disp(getenv('SLURM_JOB_NAME'))
objFuncHandle = @ODE_objective_2D_with_cebpa_runx1_PU1CEfeedback_smallestPH_0618;
ODE_fit_with_cebpa_runx1_PU1CEfeedback_PS(moreTime,objFuncHandle,nuc_destab_location_cebpa_runx1,selec_cebpa_runx1,PCn,CEPROD)
% ODE_fit_with_cebpa_runx1_PU1CEfeedback_PS(0,objFuncHandle,nuc_destab_location_cebpa_runx1,selec_cebpa_runx1); 
