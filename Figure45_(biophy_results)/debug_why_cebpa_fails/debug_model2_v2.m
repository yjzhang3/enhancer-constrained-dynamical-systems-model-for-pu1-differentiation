% how do I know if pu1 and ce needed for fitting gene expression well is
% not simulatnously good for pu1 production?

% use best-fitted pu1 and ce from model 4 to compute pu1 production rate
% assuming no runx1. 
load("best_model_ss.mat")
load('error_prop_0216.mat')
selec = 1:11;
nuc_destab_location = [0,1,0,1];
promoter_only_rescale_factor = best_pars(17);
energy_table_cebpa = energy_table([1,2,3,4,5,7,8,9,11,12,16]);

%% plug in energy pars from models to predict production rate here. Good thought, but because
% the thermodynamic model is fundamentally different, the preidcitiln will
% be different too
figure();
[final_err,D7_prob_bcell,D7_prob_ProLow,D7_prob_ProHigh] = ODE_objective_debug_model2_test_only(selec,nuc_destab_location,pu1_ss_conc,cebpa_ss_conc,energy_table_cebpa,promoter_only_rescale_factor);
bar([1,2,3],[D7_prob_bcell,D7_prob_ProLow,D7_prob_ProHigh],'FaceColor',[0.5 0.5 0.5])
hold on
errorbar([1,2,3],prodRate_ratio_all(:,2),prodRate_ratio_all(:,2)-prodRate_ratio_all(:,1),prodRate_ratio_all(:,3)-prodRate_ratio_all(:,2),'LineWidth',3,'Color','k','LineStyle','none')
legend('Predicted','Experimental')
xticks([1,2,3])
xticklabels({'B cell','pro low','pro high'})
ylabel('PU.1 production rate')
set(gca,'FontSize',20)