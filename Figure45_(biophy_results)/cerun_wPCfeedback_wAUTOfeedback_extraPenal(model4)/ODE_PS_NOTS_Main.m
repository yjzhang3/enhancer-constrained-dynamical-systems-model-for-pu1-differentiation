function ODE_PS_NOTS_Main(PCn,AUTOn,RPRC,CEPROD,moreTime)
selec_cebpa_runx1 = 1:12;
nuc_destab_location_cebpa_runx1 = [0,1,0,1,0,0];
% PCn = [2,2];
% AUTOn = [4,4];
% RPRC = [0,0];

objFuncHandle = @ODE_objective_2D_cerun_PCfeed_CEauto_smallestPH_0618;
ODE_fit_cerun_PCfeed_CEauto_PS(moreTime,objFuncHandle,nuc_destab_location_cebpa_runx1,selec_cebpa_runx1,PCn,AUTOn,RPRC,CEPROD);
