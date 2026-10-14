% plot the kueh paper proposed phenemological growth rate, converted to cell cycle length, and then
% plot the bar graph of cell cycle, and point
% out the contradiction?
clear
clc
close all
% close all
fs = 25; % overall fontsize
fs_small = 18; % legend fontsize
fa = 0.6;
lw = 3; % linewidth
ms = 13; % markersize
mk = 's'; % marker
patch_bcell_color = '#FF5800';
patch_prolow_color = '#BFBFBF';
patch_prohigh_color = 'k';
mac_color = '#2626B2';

%% plot kueh model "cell division" rate and convert to cell cycle length
pu1_totest_range_nondim = 0:0.01:1.2;
pu1_totest_range = 1:1:40000;
celldiv_rate_func = @(Cmin,Cmax,Kc,N,x) Cmin + (Cmax-Cmin)*Kc^N/(Kc^N+x^N);
Cmin = 0.02;
Cmax = 0.1;
Kc = 20000;
N = 20;
% D = 0.02;
p_ss = 37500;
norm_Kc = Kc/37500;

% there are two ways to compute this!
% rate at any dimensionless pu1 over rate at 1
cell_cycle_range_v1 = zeros(size(pu1_totest_range_nondim));
for n = 1:length(pu1_totest_range_nondim)
    cell_cycle_range_v1(n) = log(2)/celldiv_rate_func(Cmin,Cmax,norm_Kc,N,pu1_totest_range_nondim(n))/... 
    (log(2)/(celldiv_rate_func(Cmin,Cmax,norm_Kc,N,1)));
end

% rate at regular pu1 over rate at steady state
cell_cycle_range_v2 = zeros(size(pu1_totest_range));
for n = 1:length(pu1_totest_range)
    cell_cycle_range_v2(n) = log(2)/celldiv_rate_func(Cmin,Cmax,Kc,N,pu1_totest_range(n))/...
        (log(2)/(celldiv_rate_func(Cmin,Cmax,Kc,N,p_ss)));
end

%% Range of [PU1] concentration and cell cycle length as patch
load('error_prop_0216.mat');

pu1_bcell_conc_range = conc_ratio_all(1,[1,3]); % skip the median of ratio
pu1_ProLow_conc_range = conc_ratio_all(2,[1,3]);
pu1_ProHigh_conc_range = conc_ratio_all(3,[1,3]);

cellcycle_real_Bcell = cellcycle_length_ratio_all(1,[1,3]);
cellcycle_real_ProLow = cellcycle_length_ratio_all(2,[1,3]);
cellcycle_real_ProHigh = cellcycle_length_ratio_all(3,[1,3]);

f1 = [1 2 3 4];

% what if we combine production and [pu1] range?
color_patch_pu1_bcell_overall = [ pu1_bcell_conc_range(1) cellcycle_real_Bcell(1); pu1_bcell_conc_range(2) cellcycle_real_Bcell(1); pu1_bcell_conc_range(2) cellcycle_real_Bcell(2); pu1_bcell_conc_range(1) cellcycle_real_Bcell(2)]; % vertices of color patch
color_patch_pu1_ProLow_overall = [ pu1_ProLow_conc_range(1) cellcycle_real_ProLow(1); pu1_ProLow_conc_range(2) cellcycle_real_ProLow(1); pu1_ProLow_conc_range(2) cellcycle_real_ProLow(2); pu1_ProLow_conc_range(1) cellcycle_real_ProLow(2)]; % vertices of color patch
color_patch_pu1_ProHigh_overall = [ pu1_ProHigh_conc_range(1) cellcycle_real_ProHigh(1); pu1_ProHigh_conc_range(2) cellcycle_real_ProHigh(1); pu1_ProHigh_conc_range(2) cellcycle_real_ProHigh(2); pu1_ProHigh_conc_range(1) cellcycle_real_ProHigh(2)]; % vertices of color patch


%% plot experimental cell cycle normalized by macrophage removal rate
% uncertainty in ratio
f = figure('Position',[100 100 800 600]);

patch('Faces',f1,'Vertices',color_patch_pu1_bcell_overall,'FaceColor',patch_bcell_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
patch('Faces',f1,'Vertices',color_patch_pu1_ProLow_overall,'FaceColor',patch_prolow_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',patch_prohigh_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
% patch('Faces',f1,'Vertices',color_patch_pu1_mac_overall,'FaceColor','#000080','EdgeColor','none','FaceAlpha',0.75,'LineWidth',3)
hold on

plot(1,1,'LineStyle','none','Marker','pentagram','MarkerFaceColor',mac_color,'MarkerEdgeColor','none','MarkerSize',ms*1.5);

%% plot cell cycle length
axis tight
plot(pu1_totest_range_nondim,cell_cycle_range_v1,'LineWidth',lw,'Color','k','LineStyle','-')
hold on

%% overall
legend('Measured B cell','Measured pro low','Measured pro high','Macrophage (ref)','Model','Location','Northwest','FontSize',fs_small)
xlabel('\it{[pu]}')
ylabel('Cell cycle length')
xlim([0 1.2])
ylim([0 1.5])
ax = gca;
ax.Box = 'on';
set(gca,'FontSize',fs)
exportgraphics(f,'contradicting pu1-dependent cc v3.tif','Resolution',300);