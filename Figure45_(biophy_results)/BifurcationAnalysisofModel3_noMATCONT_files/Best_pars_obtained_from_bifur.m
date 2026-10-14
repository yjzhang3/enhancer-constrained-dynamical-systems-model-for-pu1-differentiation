%% 
clear
load('best_model_ss.mat','pu1_ss_conc','cebpa_ss_conc','best_pars')
x0 = [pu1_ss_conc;cebpa_ss_conc];
A = readmatrix("fixedG_1Hill.txt");
[~,I] = min(A(:,end));
p_cebp_hill = A(I,1:end-1);
p_cebp_hill_Amin = p_cebp_hill(1);
p_cebp_hill_Amax = p_cebp_hill(2);
p_cebp_hill_N = p_cebp_hill(3);
p_cebp_hill_K = p_cebp_hill(4);
% p_auto_hill_Amin = p_auto_hill(1);
% p_auto_hill_Amax = p_auto_hill(2);
% p_auto_hill_N = p_auto_hill(3);
% p_auto_hill_K = p_auto_hill(4);
MATCONTpars = [p_cebp_hill_Amin,p_cebp_hill_Amax,p_cebp_hill_N,p_cebp_hill_K,2];
ap_idx = 1;
% 
% MyContSetOptions.b = 1; 
% MyContSetOptions.maxstepsize = 0.01; % maximum step size
% MyContSetOptions.maxpt = 50; % maximum number of points
% MyContSetOptions.mintepsize = 0.00001; % min step size
% 
% odefile = @AllMATCONTFuncHandle;
% [xalpha,salpha,rangealpha,falpha] = oneparbif_weigenvalue(odefile,x0,MATCONTpars,ap_idx,MyContSetOptions);


%%
% load('CT=2_ErrVSPar_cebphill_{Amin}_Plot1.mat')
% [minerr1,I1] = min(all_err);
% disp(minerr1)
% load('CT=2_ErrVSPar_cebphill_{Amin}_Plot2.mat')
% [minerr2,I2] = min(all_err);
% disp(minerr2)
% 
% if minerr2 < minerr1
% MATCONTpars(ap_idx) = rangealpha(I2);
% else
% MATCONTpars(ap_idx) = rangealpha(I1);
% end

load('CT=2_cebphill_{Amin}_combined.mat')
[minerr,I2] = min(all_err_sorted);
MATCONTpars(ap_idx) = rangealpha_sorted(I2);

bifur_best_pars = [best_pars(1:19),MATCONTpars(1:end-1)];
objFuncHandle = @ODE_objective_2D_with_cebpa_runx1_PU1CEfeedback_smallestPH_0618;
x0_for_big_model = x0(:,1:3)';
x0_for_big_model(2,:) = [0.09,0.09];

%%
ODE_prediction_cebpa_runx1_PU1CEfeedback(objFuncHandle,1:12,bifur_best_pars,[0,1,0,1,0,0],'N',x0_for_big_model,'Model3_best_results')

