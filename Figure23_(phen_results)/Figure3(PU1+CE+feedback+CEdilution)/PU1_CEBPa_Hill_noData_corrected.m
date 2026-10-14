% this implementation is wrong!!
clear
clc
close all
load(fullfile(fileparts(pwd), 'Error propagation and Normalization', 'error_prop.mat'));
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

%% equations. now production rate = f(x)*g(x)
syms basal vmax khalf N x
syms g_basal g_vmax g_khalf g_N g_x
syms g_input_bcell g_input_prolow g_input_prohigh 
Hill_func = basal+(vmax-basal)*(x^N/(x^N+khalf^N));
g_Hill_func = g_basal+(g_vmax-g_basal)*(g_x^g_N/(g_x^g_N+g_khalf^g_N));
fg_Hill_func = Hill_func*g_Hill_func;

%%%% unknown: basal vmax khalf N g_basal g_vmax g_khalf g_N coop (8 unknowns)
% I need two more known....
% N_eval = 3;
% basal_eval = 0.001;
% fg_Hill_func = subs(fg_Hill_func,[N,basal],[N_eval,basal_eval]);

%%%% since you plug in values of x to get a new equation, x is no longer an
% unknown we want to solve!!
fg_Hill_func_Bcell = subs(fg_Hill_func,[x,g_x],[conc_ratio_all(1,2),g_input_bcell]);
fg_Hill_func_ProLow = subs(fg_Hill_func,[x,g_x],[conc_ratio_all(2,2),g_input_prolow]);
fg_Hill_func_ProHigh = subs(fg_Hill_func,[x,g_x],[conc_ratio_all(3,2),g_input_prohigh]);
fg_Hill_func_mac = subs(fg_Hill_func,[x,g_x],[1,1]);

% now seven equations
eqn1 = fg_Hill_func_Bcell == remRate_constant_ratio_all(1,2)*conc_ratio_all(1,2);
eqn2 = fg_Hill_func_ProLow == remRate_constant_ratio_all(2,2)*conc_ratio_all(2,2);
eqn3 = fg_Hill_func_ProHigh == remRate_constant_ratio_all(3,2)*conc_ratio_all(3,2);
eqn4 = fg_Hill_func_mac == 1*1;
eqn5 = fg_Hill_func_Bcell == prodRate_ratio_all(1,2);
eqn6 = fg_Hill_func_ProLow == prodRate_ratio_all(2,2);
eqn7 = fg_Hill_func_ProHigh == prodRate_ratio_all(3,1);

%% solve
% assume(basal >0)
% % assume(n, 'integer')
% assume([x basal], 'real')       % Coefficients are real numbers  
% S = solve([eqn1], [x],'ReturnConditions',true);
S = vpasolve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7], [basal,vmax,khalf,N,g_basal,g_vmax,g_khalf,g_N,g_input_bcell,g_input_prolow,g_input_prohigh]);
% number of unknowns > number of equations

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

subplot(2,2,1)
% uncertainty in ratio
patch('Faces',f1,'Vertices',color_patch_pu1_bcell_overall,'FaceColor',patch_bcell_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
plot(conc_ratio_all(1,2),prodRate_ratio_all(1,2),'LineStyle','none','Marker','+','MarkerFaceColor',patch_bcell_color,'MarkerEdgeColor','#FF5800','MarkerSize',ms*0.8,'LineWidth',6);

subplot(2,2,2)
patch('Faces',f1,'Vertices',color_patch_pu1_ProLow_overall,'FaceColor',patch_prolow_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
plot(conc_ratio_all(2,2),prodRate_ratio_all(2,2),'LineStyle','none','Marker','+','MarkerFaceColor',patch_prolow_color,'MarkerEdgeColor','#BFBFBF','MarkerSize',ms*0.8,'LineWidth',6);

subplot(2,2,3)
patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',patch_prohigh_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
plot(conc_ratio_all(3,2),prodRate_ratio_all(3,2),'LineStyle','none','Marker','+','MarkerFaceColor',patch_prohigh_color,'MarkerEdgeColor','k','MarkerSize',ms*0.8,'LineWidth',6);

subplot(2,2,4)
plot(1,1,'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor',remv_color{4},'MarkerSize',ms)


pu1_totest_conc = 0:0.1:1.2;
for ct = 1:3
    % same color, different mk
    % plot(pu1_totest_conc,remRate_constant_ratio_all(ct,2)*pu1_totest_conc,'LineWidth',lw,'Color',remv_color,'Marker',remv_celltype_marker{ct},'MarkerEdgeColor',remv_color,'MarkerFaceColor','none','MarkerSize',ms)
    subplot(2,2,ct)
    % same mk, different color
    plot(pu1_totest_conc,remRate_constant_ratio_all(ct,2)*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{ct})

    hold on
end

% plot(pu1_totest_conc,1*pu1_totest_conc,'LineWidth',lw*0.8,'Color',remv_color,'Marker',remv_celltype_marker{4},'MarkerEdgeColor',remv_color,'MarkerFaceColor','none','MarkerSize',ms)
subplot(2,2,4)
plot(pu1_totest_conc,1*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{4})



ax = findall(gcf, 'Type', 'axes');
xlim(ax,[0 1.2])
ylim(ax,[0 3])
set(ax, 'FontSize', fs)
% set(ax(4),'YScale','log','XScale','log')
















