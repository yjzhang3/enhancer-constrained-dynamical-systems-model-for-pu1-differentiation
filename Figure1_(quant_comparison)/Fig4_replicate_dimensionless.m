clear
clc
close all
fa = 0.6; % facealpha for patch
fs = 25; % fontsize
lg_fs = 18;
ms = 13; % marker size
lw = 3;
patch_bcell_color = '#FF5800';
patch_prolow_color = '#BFBFBF';
patch_prohigh_color = 'k';
kueh_pro_color = 'c';
mac_color = '#2626B2';
ss_ylim = 1.2;
prod_ylim = 2.5;

%% Range of [PU1] concentration in which we find intersection (color patch)
load('error_prop_0216.mat');

pu1_bcell_conc_range = conc_ratio_all(1,[1,3]); % skip the median of ratio
pu1_ProLow_conc_range = conc_ratio_all(2,[1,3]);
pu1_ProHigh_conc_range = conc_ratio_all(3,[1,3]);

pu1_bcell_prod_rate_range = prodRate_ratio_all(1,[1,3]);
pu1_ProLow_prod_rate_range = prodRate_ratio_all(2,[1,3]);
pu1_ProHigh_prod_rate_range = prodRate_ratio_all(3,[1,3]);

f1 = [1 2 3 4];

% what if we combine production and [pu1] range?
color_patch_pu1_bcell_overall = [ pu1_bcell_conc_range(1) pu1_bcell_prod_rate_range(1); pu1_bcell_conc_range(2) pu1_bcell_prod_rate_range(1); pu1_bcell_conc_range(2) pu1_bcell_prod_rate_range(2); pu1_bcell_conc_range(1) pu1_bcell_prod_rate_range(2)]; % vertices of color patch
color_patch_pu1_ProLow_overall = [ pu1_ProLow_conc_range(1) pu1_ProLow_prod_rate_range(1); pu1_ProLow_conc_range(2) pu1_ProLow_prod_rate_range(1); pu1_ProLow_conc_range(2) pu1_ProLow_prod_rate_range(2); pu1_ProLow_conc_range(1) pu1_ProLow_prod_rate_range(2)]; % vertices of color patch
color_patch_pu1_ProHigh_overall = [ pu1_ProHigh_conc_range(1) pu1_ProHigh_prod_rate_range(1); pu1_ProHigh_conc_range(2) pu1_ProHigh_prod_rate_range(1); pu1_ProHigh_conc_range(2) pu1_ProHigh_prod_rate_range(2); pu1_ProHigh_conc_range(1) pu1_ProHigh_prod_rate_range(2)]; % vertices of color patch

%% replicate figure 4 in paper
A_max = 1500;
A_min = 165;
C_max = 0.1;
C_min = 0.02;
M = 5;
N = 20;
p_ss  = 37500;
Kp = 6000; % it doesn't matter which steady state you use for normalization, there are always 3 stable ss (one of them would have a value of 1, which 
% is the one you chose for normalization)
Kc = 20000;
D = 0.02;
pu1_totest_range = 0:1/40000:ss_ylim; % this is dimensionless concentration;
% adjust this upper bound accordingly if the highest stable steady state is not used for normalization
Ap_all = zeros(length(pu1_totest_range),1);
CDp_all = zeros(length(pu1_totest_range),1);
for yt = 1:length(pu1_totest_range)
%     Ap_all(yt) = production_rate(p(yt),A_min,A_max,M,Kp/p_ss)/production_rate(p_ss,A_min,A_max,M,Kp);
%     CDp_all(yt) = (dilution_rate(p(yt),C_max,C_min,N,Kc/p_ss)+D)/(dilution_rate(p_ss,C_max,C_min,N,Kc)+D)*p(yt);
    % same thing
%     Ap_all(yt) = production_rate(p(yt)*p_ss,A_min,A_max,M,Kp)/production_rate(p_ss,A_min,A_max,M,Kp);
%     CDp_all(yt) = (dilution_rate(p(yt)*p_ss,C_max,C_min,N,Kc)+D)/(dilution_rate(p_ss,C_max,C_min,N,Kc)+D)*p(yt);
    % same thing
    Ap_all(yt) = production_rate(pu1_totest_range(yt)*p_ss,A_min,A_max,M,Kp)/production_rate(1,A_min,A_max,M,Kp/p_ss);
    CDp_all(yt) = (dilution_rate(pu1_totest_range(yt)*p_ss,C_max,C_min,N,Kc)+D)/(dilution_rate(1,C_max,C_min,N,Kc/p_ss)+D)*pu1_totest_range(yt);
end

ind = find(abs(Ap_all-CDp_all) < 0.0001); % manually find intersections
concise_ind_arr = extract_intersection_index(ind,40);
%% Figure 1A
% close all
f = figure('Position',[100 100 800 600]);
axis tight

% uncertainty in ratio
% patch('Faces',f1,'Vertices',color_patch_pu1_bcell_overall,'FaceColor',patch_bcell_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
% hold on
% patch('Faces',f1,'Vertices',color_patch_pu1_ProLow_overall,'FaceColor',patch_prolow_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
% hold on
% patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',patch_prohigh_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
% hold on
% % patch('Faces',f1,'Vertices',color_patch_pu1_mac_overall,'FaceColor','#000080','EdgeColor','none','FaceAlpha',0.75,'LineWidth',3)
% hold on

% median in ratio
% plot(conc_ratio_all(1,2),prodR_ratio_all(1,2),'LineStyle','none','Marker','+','MarkerFaceColor','#FF5800','MarkerEdgeColor','#FF5800','MarkerSize',ms*0.8,'LineWidth',6);
% hold on
% plot(conc_ratio_all(2,2),prodR_ratio_all(2,2),'LineStyle','none','Marker','+','MarkerFaceColor','#BFBFBF','MarkerEdgeColor','#BFBFBF','MarkerSize',ms*0.8,'LineWidth',6);
% hold on
% plot(conc_ratio_all(3,2),prodR_ratio_all(3,2),'LineStyle','none','Marker','+','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',ms*0.8,'LineWidth',6);

plot(pu1_totest_range,Ap_all,'LineWidth',lw,'Color',[0.00 0.75 0.65])
hold on
plot(pu1_totest_range,CDp_all,'LineWidth',lw,'Color',[0.70 0.40 0.85])
ylim([0 prod_ylim])
xlim([0 pu1_totest_range(end)])

% xline(p(1475),'LineStyle',':','LineWidth',3,'Color','#FF5800')
% hold on
% xline(p(13002),'LineStyle',':','LineWidth',3,'Color','#808080')
% hold on
% xline(p(40000),'LineStyle',':','LineWidth',3,'Color','#000080')

plot(pu1_totest_range(concise_ind_arr(1)),Ap_all(concise_ind_arr(1)),'Marker','o','MarkerSize',ms,'LineStyle','none','MarkerEdgeColor',patch_bcell_color,'MarkerFaceColor',patch_bcell_color);
plot(pu1_totest_range(concise_ind_arr(3)),Ap_all(concise_ind_arr(3)),'Marker','o','MarkerSize',ms,'LineStyle','none','MarkerEdgeColor',kueh_pro_color,'MarkerFaceColor',kueh_pro_color);
plot(pu1_totest_range(concise_ind_arr(5)),Ap_all(concise_ind_arr(5)),'Marker','pentagram','MarkerSize',ms*1.5,'LineStyle','none','MarkerEdgeColor',mac_color,'MarkerFaceColor',mac_color);

% legend('Data-constrained B cell','Data-constrained pro low','Data-constrained pro high','Model-predicted prod. rate','Model-predicted remv. rate','Location','northeast','fontsize',lg_fs)
legend('Production rate','Removal rate','B cell','Progenitor cell','Macrophage (ref)','Location','northeast','fontsize',lg_fs)


% xlabel('Dimensionless Protein Concentration \it{[pu]}')
xlabel('\it{[pu]}')
ylabel('Rate')

ax = gca;
ax.Box = 'on';

set(gca,'FontSize',fs)
exportgraphics(f,'replicate_fig4_dimensionless_v2.tif','Resolution',300);
% saveas(gcf,'replicate_fig4_dimensionless.tif')

%% also bring up our results presentation for the first time
% show_results_real_vs_pred(p(concise_ind_arr([1,3,3])),...
%     Ap_all(concise_ind_arr([1,3,3])),'Kueh_paper_')
% 
% %%
% show_results_real_vs_pred_v2(p(concise_ind_arr([1,3,3])),...
%     Ap_all(concise_ind_arr([1,3,3])),'Kueh_paper_')
% close all

%% Figure 1B
show_results_real_vs_pred_v3_Kueh(pu1_totest_range(concise_ind_arr([1,3,3])),Ap_all(concise_ind_arr([1,3,3])),'Kueh_paper')