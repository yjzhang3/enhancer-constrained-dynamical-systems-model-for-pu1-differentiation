% plug less rigorous method-yielding results into ODE
%% read best pars from less rigious method
clear
clc
selec = 1:12;
nuc_destab_location = [0,1,0,1,0,0]; % there is no RUNX1 in the model
lw = 3; % linewidth
ms = 12; % marker size
fs = 25; % fontsize
lg_fs = 19;
fa = 0.6; % facealpha of constraint patch (not growth rate patch)
% color for constraint patch
bcellcolor = '#FF5800';
prolowcolor = '#BFBFBF';
prohighcolor = 'k';
maccolor = '#2626B2';

final_fn = "local_fit_with_cebpa_runx1_general_PU1CEfeedback_0401_PattrC_2N_maxV1_123456789101112_010100.txt";
A = readmatrix(final_fn);
arri = A(:,end);
sorted_arri = sort(arri);
minerr = sorted_arri(1); % pick the smallest or second smallest?
minind = find(arri==minerr);
best_pars = A(minind,1:end-1);


[final_er,sim_d_final,real_d_final,energy_table,...
    concentration_table,concentration_table_Bcell,...
    concentration_table_ProLow,concentration_table_ProHigh,all_sublist_str,...
    norm_D7_prob_bcell,norm_D7_prob_ProLow,norm_D7_prob_ProHigh,...
    p_cebp_hill,pu1_prod_slope_celltype] = objective_with_cebpa_runx1_general_PU1CEfeedback_0401([1:12],best_pars,0,0,0,nuc_destab_location);

model_predicted_Bcell_pu1 = concentration_table_Bcell(2);
model_predicted_ProLow_pu1 = concentration_table_ProLow(2);
model_predicted_ProHigh_pu1 = concentration_table_ProHigh(2);
model_predicted_Bcell_cebpa = concentration_table_Bcell(3);
model_predicted_ProLow_cebpa = concentration_table_ProLow(3);
model_predicted_ProHigh_cebpa = concentration_table_ProHigh(3);

pu1_ss = [model_predicted_Bcell_pu1,model_predicted_ProLow_pu1,model_predicted_ProHigh_pu1];
ce_ss = [model_predicted_Bcell_cebpa,model_predicted_ProLow_cebpa,model_predicted_ProHigh_cebpa];

%% use predicted pu1 and cebpa conc to simulate ODE
clc
p = best_pars(1,[1:17,21:26]);
pu1_x0 = [model_predicted_Bcell_pu1,model_predicted_ProLow_pu1,model_predicted_ProHigh_pu1];
pu1_x0(2) = 0.15;
cebpa_x0 = [model_predicted_Bcell_cebpa,model_predicted_ProLow_cebpa,model_predicted_ProHigh_cebpa];
% cebpa_x0(2) = 0.01;

x0 = [pu1_x0',cebpa_x0'];



[final_err,pu1_ss_conc,cebpa_ss_conc,runx1_ss_conc,pu1_ss_prod,...
    real_d_final,sim_d_final,D7_valid_sublist,energy_table,...
    Bcell_over_time2,ProLow_over_time2,ProHigh_over_time2,tB2,tPL2,tPH2] = ...
    ODE_objective_2D_with_cebpa_runx1_PU1CEfeedback_smallestPH_0401(selec,p,nuc_destab_location,x0);

% the original steady state leads to a good error, but steady state of
% prolow is way too high. ODE solver with the wrong initial condition cna't
% improve. So I'd rather lower the initial condition and start with a much
% higher error
%% plot time dynamics
f = figure('Position',[100 100 800 600]);
axis tight
plot(tB2,Bcell_over_time2(:,1),'LineWidth',lw,'Color',bcellcolor)
hold on
plot(tB2,Bcell_over_time2(:,2),'LineWidth',lw,'Color',bcellcolor,'LineStyle','--')
hold on
plot(tPL2,ProLow_over_time2(:,1),'LineWidth',lw,'Color',prolowcolor)
hold on
plot(tPL2,ProLow_over_time2(:,2),'LineWidth',lw,'Color',prolowcolor,'LineStyle','--')
hold on
plot(tPH2,ProHigh_over_time2(:,1),'LineWidth',lw,'Color',prohighcolor)
hold on
plot(tPH2,ProHigh_over_time2(:,2),'LineWidth',lw,'Color',prohighcolor,'LineStyle','--')

xlabel('t')
ylabel('\it{[pu]}')
ylim([0 1])
legend('B cell','','pro low','','pro high','','FontSize',lg_fs)
set(gca,'FontSize',fs)

%% determine whether to save this x0
save('wRUNX1_noODEsim_ce_conc_pars_0401.mat',"pu1_ss","ce_ss","x0")