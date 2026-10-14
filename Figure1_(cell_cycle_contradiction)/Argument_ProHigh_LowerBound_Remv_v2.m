%% this file is to show the reason why we decided to use the lower bound of pro
% high removal rate for optimization for all models
clear
clc
close all
load('error_prop_0216.mat');

lw = 3;
ms = 15;
fa = 0.6;
fs = 25;
fs_small = 16;

PH_removcolor = '#808080';
PH_patchcolor = 'k';
ss_ylim = 1.2;
prod_ylim = 2.5;

pu1_totest_conc = 0:0.0001:ss_ylim;

%% Range of [PU1] concentration in which we find intersection (color patch)
pu1_ProHigh_conc_range = conc_ratio_all(3,[1,3]);
pu1_ProHigh_prod_rate_range = prodRate_ratio_all(3,[1,3]);

f1 = [1 2 3 4];

% what if we combine production and [pu1] range?
color_patch_pu1_ProHigh_overall = [ pu1_ProHigh_conc_range(1) pu1_ProHigh_prod_rate_range(1); pu1_ProHigh_conc_range(2) pu1_ProHigh_prod_rate_range(1); pu1_ProHigh_conc_range(2) pu1_ProHigh_prod_rate_range(2); pu1_ProHigh_conc_range(1) pu1_ProHigh_prod_rate_range(2)]; % vertices of color patch

%% plot
close all
f = figure('Position',[293,212,1571,709]);
axis tight

subplot(1,2,1)
patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',PH_patchcolor,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
% 
% % plot_range_growthrate_patch(pu1_totest_conc,remRate_constant_ratio_all(3,3)*pu1_totest_conc,remRate_constant_ratio_all(3,1)*pu1_totest_conc,hex2rgb(PH_removcolor))
% plot_range_growthrate_patch(pu1_totest_conc,remRate_constant_ratio_all(3,3)*pu1_totest_conc,remRate_constant_ratio_all(3,1)*pu1_totest_conc,[1 0 0])
% hold on
plot(pu1_totest_conc,remRate_constant_ratio_all(3,2)*pu1_totest_conc,'LineWidth',lw,'Color','r')
% hold on
% plot(pu1_totest_conc,remRate_constant_ratio_all(3,1)*pu1_totest_conc,'LineWidth',lw,'Color','b')
% hold on

ind_original = find(abs(pu1_totest_conc- conc_ratio_all(3,1)) < 1e-4); % because of the resolution of pu1_totest_conc,
% this just helps finds the concentration that is the lower bound of PU1gfp

% plot where the median removal rate will be if constrained by pu1 conc of
% pro high
plot(pu1_ProHigh_conc_range(1),remRate_constant_ratio_all(3,2)*pu1_ProHigh_conc_range(1),'LineStyle','none','Marker','^','LineWidth',lw,'MarkerFaceColor','none','MarkerEdgeColor','r','MarkerSize',ms)
hold on

% plot where the lower bound removal rate will be if constrained by pu1 conc of
% pro high
% yline(remRate_constant_ratio_all(3,1)*pu1_totest_conc(ind_original(1)),'Color','b','LineStyle','--','LineWidth',lw)

% final touch
xlim([0 ss_ylim])
ylim([0 prod_ylim])
xlabel('\it{[pu]}')
ylabel('PU.1 synthesis rate')

legend('Real Pro High','Removal rate constrained by median cc length',sprintf('Min. synthesis rate at SS = %.3f',remRate_constant_ratio_all(3,2)*pu1_ProHigh_conc_range(1)), ...
    'Location','northoutside','fontsize',fs_small)

ax = gca;
ax.Box = 'on';
ax.FontSize = fs;
%%
subplot(1,2,2)
patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',PH_patchcolor,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on

plot(pu1_totest_conc,remRate_constant_ratio_all(3,1)*pu1_totest_conc,'LineWidth',lw,'Color','b')
hold on

% plot where the lower bound removal rate will be if constrained by pu1 conc of
% pro high
rmv_func = remRate_constant_ratio_all(3,1)*pu1_totest_conc;
mask1 = rmv_func<1;
mask2 = pu1_totest_conc>pu1_totest_conc(ind_original(1));

plot(pu1_totest_conc(logical(mask1.*mask2)), ...
    rmv_func(logical(mask1.*mask2)),'LineStyle','none','MarkerFaceColor','none','MarkerEdgeColor','b','LineWidth',lw,'Marker','^','MarkerSize',ms)

legend('Real Pro High','Removal rate constrained by 25% quantile cc length',...
    'All possible SS with synthesis rate < 1',...
    'Location','northoutside','fontsize',fs_small)

xlim([0 ss_ylim])
ylim([0 prod_ylim])
xlabel('\it{[pu]}')
ylabel('PU.1 synthesis rate')


ax = gca;
ax.Box = 'on';
ax.FontSize = fs;

exportgraphics(f,'ProHigh_LowerBound_Remv.tif','Resolution',300);