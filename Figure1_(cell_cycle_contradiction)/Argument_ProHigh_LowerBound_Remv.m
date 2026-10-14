%% this file is to show the reason why we decided to use the lower bound of pro
% high removal rate for optimization for all models
clear
clc
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
f = figure('Position',[100 100 800 600]);
axis tight

patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',PH_patchcolor,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on

% plot_range_growthrate_patch(pu1_totest_conc,remRate_constant_ratio_all(3,3)*pu1_totest_conc,remRate_constant_ratio_all(3,1)*pu1_totest_conc,hex2rgb(PH_removcolor))
plot_range_growthrate_patch(pu1_totest_conc,remRate_constant_ratio_all(3,3)*pu1_totest_conc,remRate_constant_ratio_all(3,1)*pu1_totest_conc,[1 0 0])
hold on
plot(pu1_totest_conc,remRate_constant_ratio_all(3,2)*pu1_totest_conc,'LineWidth',lw,'Color','r')
hold on
plot(pu1_totest_conc,remRate_constant_ratio_all(3,1)*pu1_totest_conc,'LineWidth',lw,'Color','b')
hold on

ind_original = find(abs(pu1_totest_conc- conc_ratio_all(3,1)) < 1e-4); % because of the resolution of pu1_totest_conc,
% this just helps finds the concentration that is the lower bound of PU1gfp

% plot where the median removal rate will be if constrained by pu1 conc of
% pro high
yline(remRate_constant_ratio_all(3,2)*pu1_totest_conc(ind_original(1)),'Color','r','LineStyle','--','LineWidth',lw)
hold on

% plot where the lower bound removal rate will be if constrained by pu1 conc of
% pro high
yline(remRate_constant_ratio_all(3,1)*pu1_totest_conc(ind_original(1)),'Color','b','LineStyle','--','LineWidth',lw)

%% final touch
xlim([0 ss_ylim])
ylim([0 prod_ylim])
xlabel('\it{[pu]}')
ylabel('PU.1 synthesis rate')

legend('Real Pro High','','Median removal rate','25% quantile removal rate','SS synthesis rate given median removal rate', ...
    'SS synthesis rate required given 25% quantile removal rate',...
    'Location','northoutside','fontsize',fs_small)

ax = gca;
ax.Box = 'on';
ax.FontSize = fs;

exportgraphics(f,'ProHigh_LowerBound_Remv.tif','Resolution',300);