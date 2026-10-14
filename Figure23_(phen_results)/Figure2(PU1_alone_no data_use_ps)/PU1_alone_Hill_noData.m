clear
clc
close all
load(fullfile(fileparts(pwd), 'Figure1A', 'error_prop.mat'));
% remv_color = [0.70 0.40 0.85];
remv_color = {'#FF5800','#BFBFBF','#808080','#2626B2'}; % or different color for different ct
% remv_celltype_marker = {'+','d','*','s'}; % either different marker for different ct
lw = 3;
ms = 10;
fa = 0.6;
fs = 19;
patch_bcell_color = '#FF5800';
patch_prolow_color = '#BFBFBF';
patch_prohigh_color = 'k';

%% Range of [PU1] concentration in which we find intersection (color patch)
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

f = figure('Position',[100 100 800 600]);
axis tight

% uncertainty in ratio
patch('Faces',f1,'Vertices',color_patch_pu1_bcell_overall,'FaceColor',patch_bcell_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
patch('Faces',f1,'Vertices',color_patch_pu1_ProLow_overall,'FaceColor',patch_prolow_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',patch_prohigh_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
% patch('Faces',f1,'Vertices',color_patch_pu1_mac_overall,'FaceColor','#000080','EdgeColor','none','FaceAlpha',0.75,'LineWidth',3)
hold on
plot(1,1,'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor',remv_color{4},'MarkerSize',ms)

% median in ratio
plot(conc_ratio_all(1,2),prodRate_ratio_all(1,2),'LineStyle','none','Marker','+','MarkerFaceColor',patch_bcell_color,'MarkerEdgeColor','#FF5800','MarkerSize',ms*0.8,'LineWidth',6);
hold on
plot(conc_ratio_all(2,2),prodRate_ratio_all(2,2),'LineStyle','none','Marker','+','MarkerFaceColor',patch_prolow_color,'MarkerEdgeColor','#BFBFBF','MarkerSize',ms*0.8,'LineWidth',6);
hold on
plot(conc_ratio_all(3,2),prodRate_ratio_all(3,2),'LineStyle','none','Marker','+','MarkerFaceColor',patch_prohigh_color,'MarkerEdgeColor','k','MarkerSize',ms*0.8,'LineWidth',6);

xlim([0 1.2])
ylim([0 2])
set(gca,'FontSize',fs)

%% equations
syms basal vmax khalf N x

Hill_func = basal+(vmax-basal)*(x^N/(x^N+khalf^N));
Hill_func_Bcell = subs(Hill_func,[x],[conc_ratio_all(1,2)]);
Hill_func_ProLow = subs(Hill_func,[x],[conc_ratio_all(2,2)]);
Hill_func_ProHigh = subs(Hill_func,[x],[conc_ratio_all(3,2)]);
Hill_func_mac = subs(Hill_func,[x],[1]);

% four equatinos can either be let prod = rmv
% doesn't work because 1) prohigh rmv rate > prolow rmv rate, so prod
% intersects prohigh first and then prolow; 2) however, pu1 conc at prohigh
% must be larger than conc at prolow (i.e., intersections must happen for
% prolow first)
% eqn1 = Hill_func_Bcell == remRate_constant_ratio_all(1,2)*conc_ratio_all(1,2);
% eqn2 = Hill_func_ProLow == remRate_constant_ratio_all(2,2)*conc_ratio_all(2,2);
% eqn3 = Hill_func_ProHigh == remRate_constant_ratio_all(3,2)*conc_ratio_all(3,2);
% eqn4 = Hill_func_mac == 1*1;

% four equatinos can be let prod = data constrained prod
% doesn't work if using median of production rate because prohigh prod rate
% is higher than mac while prohigh conc is lower than mac (hill is
% monotonic). Unless we use lower bound of prohigh
eqn1 = Hill_func_Bcell == prodRate_ratio_all(1,2);
eqn2 = Hill_func_ProLow == prodRate_ratio_all(2,2);
eqn3 = Hill_func_ProHigh == prodRate_ratio_all(3,2);
eqn4 = Hill_func_mac == 1*1;

%% solve
assume(x,'real')
assume(x,'positive')
% assume(basal >0)
% % assume(n, 'integer')
% assume([x basal], 'real')       % Coefficients are real numbers  
% S = solve([eqn1], [x],'ReturnConditions',true);
S = vpasolve([eqn1,eqn2,eqn3,eqn4], [basal,vmax,khalf,N]);


%% why no solution? plot it first
f = figure('Position',[100 100 800 600]);
axis tight
pu1_totest_conc = 0:0.1:1.2;
for ct = 1:3
    % same color, different mk
    % plot(pu1_totest_conc,remRate_constant_ratio_all(ct,2)*pu1_totest_conc,'LineWidth',lw,'Color',remv_color,'Marker',remv_celltype_marker{ct},'MarkerEdgeColor',remv_color,'MarkerFaceColor','none','MarkerSize',ms)
    
    % same mk, different color
    plot(pu1_totest_conc,remRate_constant_ratio_all(ct,2)*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{ct})

    hold on
end

% plot(pu1_totest_conc,1*pu1_totest_conc,'LineWidth',lw*0.8,'Color',remv_color,'Marker',remv_celltype_marker{4},'MarkerEdgeColor',remv_color,'MarkerFaceColor','none','MarkerSize',ms)

plot(pu1_totest_conc,1*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{4})

xlim([0 1.2])
ylim([0 3])
set(gca,'FontSize',fs)
