function ODE_PS_NOTS_Main(moreTime,PCn,CEPROD)
selec = 1:11;
nuc_destab_location = [0,1,0,1];

%objFuncHandle = @objective_with_cebpa_runx1_general_PU1CEfeedback_0401;
%fit_with_cebpa_runx1_general_PU1CEfeedback(objFuncHandle,nuc_destab_location_cebpa_runx1,selec_cebpa_runx1); 
% disp(getenv('SLURM_JOB_NAME'))
objFuncHandle = @ODE_objective_2D_cebpa_PCnoAUTOfeedback_smallestPHdilute_0618;
ODE_fit_PS_cebpa_PCAUTOfeedback(moreTime,objFuncHandle,nuc_destab_location,selec,PCn,CEPROD);
