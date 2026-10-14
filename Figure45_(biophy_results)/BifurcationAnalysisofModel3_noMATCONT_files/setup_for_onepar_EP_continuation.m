function [rangealpha,all_err] = setup_for_onepar_EP_continuation(ct,odefile,MyContSetOptions,...
    search_dir,ap_idx,MATCONTpars,PlotOptions,plot_color_idx,ax1,ax2,OutDirName)
% search_dir = 0 or 1 (0 = forward, backward = 1)
% ct = cell type
% ap_idx = index of active par in pars array
% ap_initial_par_val = the initial value of active par

load('color_pal.mat','cols');
cols = cols([2,3,4,8],:);

selec = 1:12;
nuc_destab_location = [0,1,0,1,0,0];

%% set up equilibrium point
MyContSetOptions.b = search_dir; 
% MATCONTpars(ap_idx) = ap_initial_par_val;
load('best_model_ss.mat','pu1_ss_conc','cebpa_ss_conc','best_pars')
x0 = [pu1_ss_conc(ct);cebpa_ss_conc(ct)];

%% equilibrium continuation
[xalpha,salpha,rangealpha,falpha] = oneparbif_weigenvalue(odefile,x0,MATCONTpars,ap_idx,MyContSetOptions);


%% plot bifrucation diagram for pu1
ivar = 1; % pu1 is the first state variable
plotonepar_wptName(xalpha,salpha,rangealpha,falpha,ivar,PlotOptions,cols(plot_color_idx,:),ax1);

%% plot bifrucation diagram for cebpa
ivar = 2; % cebpa is the first state variable
plotonepar_wptName(xalpha,salpha,rangealpha,falpha,ivar,PlotOptions,cols(plot_color_idx,:),ax2);

%% based on the explored parameter, adjust the rest of the hill pars to ensure best fit;
% then, test error of my model for all found equilibrium 
if exist('xalpha', 'var') == 1 
x0_for_big_model = [pu1_ss_conc;cebpa_ss_conc]';
all_err = zeros(size(rangealpha));
for oo = 1:length(rangealpha)
    % adjusted_best_pars_oneHill = online_adjust_HillPars(ap_idx,rangealpha(oo));
    MATCONTpars(ap_idx) = rangealpha(oo);
    full_pars = [best_pars(1:19),MATCONTpars];
    
    % update initial condition
    x0_for_big_model(ct,:) = xalpha(1:2,oo);

    all_err(oo) = ...
    ODE_objective_2D_with_cebpa_runx1_PU1CEfeedback_smallestPH_0618(selec,full_pars,nuc_destab_location,'N',x0_for_big_model);
end
end

%% find minerr, and predict on the fly
[minerr,I] = min(all_err);
disp(minerr)
MATCONTpars(ap_idx) = rangealpha(I);
bifur_best_pars = [best_pars(1:19),MATCONTpars];
objFuncHandle = @ODE_objective_2D_with_cebpa_runx1_PU1CEfeedback_smallestPH_0618;
x0_for_big_model(ct,:) = xalpha(1:2,I);
ODE_prediction_cebpa_runx1_PU1CEfeedback(objFuncHandle,1:12,bifur_best_pars,[0,1,0,1,0,0],'N',x0_for_big_model,OutDirName)
