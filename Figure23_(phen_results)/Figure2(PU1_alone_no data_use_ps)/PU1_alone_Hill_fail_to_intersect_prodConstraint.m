% four equatinos can be let prod = data constrained prod
% doesn't work if using median of production rate because prohigh prod rate
% is higher than mac while prohigh conc is lower than mac (hill is
% monotonic). Unless we use lower bound of prohigh

%%
clear
clc
close all
load(fullfile(fileparts(pwd), 'Error propagation and Normalization', 'error_prop.mat'));
prod_color = [0.00 0.75 0.65];
% remv_color = [0.70 0.40 0.85];
remv_color = {'#FF5800','#BFBFBF','#808080','#2626B2'}; % or different color for different ct
% remv_celltype_marker = {'+','d','*','s'}; % either different marker for different ct
lw = 3;
ms = 15;
fa = 0.6;
fs = 25;
fs_small = 16;
patch_bcell_color = '#FF5800';
patch_prolow_color = '#BFBFBF';
patch_prohigh_color = 'k';
ss_ylim = 1.2;
prod_ylim = 2.5;

%% equations
syms basal vmax khalf N x

Hill_func = basal+(vmax-basal)*(x^N/(x^N+khalf^N));
Hill_func_Bcell = subs(Hill_func,[x],[conc_ratio_all(1,2)]);
Hill_func_ProLow = subs(Hill_func,[x],[conc_ratio_all(2,2)]);
Hill_func_ProHigh = subs(Hill_func,[x],[conc_ratio_all(3,2)]);
Hill_func_mac = subs(Hill_func,[x],[1]);

eqn1 = Hill_func_Bcell == prodRate_ratio_all(1,2);
eqn2 = Hill_func_ProLow == prodRate_ratio_all(2,2);
eqn3 = Hill_func_ProHigh == prodRate_ratio_all(3,1); % let's change this to lower bound of constraint to show a solution cna be found only if 
% constrained production rate is monotonic as well
eqn4 = Hill_func_mac == 1*1;

%% solve
assume(x,'real')
assume(x,'positive')
% assume(basal >0)
% % assume(n, 'integer')
% assume([x basal], 'real')       % Coefficients are real numbers  
% S = solve([eqn1], [x],'ReturnConditions',true);
S = vpasolve([eqn1,eqn2,eqn3,eqn4], [basal,vmax,khalf,N]);

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

% median in ratio
% plot(conc_ratio_all(1,2),prodRate_ratio_all(1,2),'LineStyle','none','Marker','+','MarkerFaceColor',patch_bcell_color,'MarkerEdgeColor','#FF5800','MarkerSize',ms*0.8,'LineWidth',6);
% hold on
% plot(conc_ratio_all(2,2),prodRate_ratio_all(2,2),'LineStyle','none','Marker','+','MarkerFaceColor',patch_prolow_color,'MarkerEdgeColor','#BFBFBF','MarkerSize',ms*0.8,'LineWidth',6);
% hold on
% plot(conc_ratio_all(3,2),prodRate_ratio_all(3,2),'LineStyle','none','Marker','+','MarkerFaceColor',patch_prohigh_color,'MarkerEdgeColor','k','MarkerSize',ms*0.8,'LineWidth',6);
% hold on
plot(1,1,'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',remv_color{4},'MarkerSize',ms)

%% overlay monotonic production rate found by sovling four equations
basal_sol = eval(S(1).basal);
vmax_sol =  eval(S(1).vmax);
khalf_sol =  eval(S(1).khalf);
N_sol = eval(S(1).N);
Hill_pars = [basal_sol,vmax_sol,N_sol,khalf_sol];
save("prodConstraint_hillpars.mat",'Hill_pars');

pu1_totest_conc = 0:0.01:ss_ylim;
anyn_Hill = @(basal,vmax,khalf,N,x) basal+(vmax-basal)*(x^N/(x^N+khalf^N));
prod_totest_conc = zeros(size(pu1_totest_conc));
for i = 1:length(pu1_totest_conc)
    prod_totest_conc(i) = anyn_Hill(basal_sol,vmax_sol,khalf_sol,N_sol,pu1_totest_conc(i));
end

plot(pu1_totest_conc,prod_totest_conc,'LineWidth',lw,'Color',prod_color)

xlim([0 ss_ylim])
ylim([0 prod_ylim])

xlabel('\it{[pu]}')
ylabel('Normalized rate')
legend('Data-constrained B cell','Data-constrained pro low','Data-constrained pro high','Macrophage','Phen. prod. rate','Location','northeast','fontsize',fs_small)

ax = gca;
ax.Box = 'on';
set(gca,'FontSize',fs)
exportgraphics(f,'PU1_alone_Hill_fail_to_intersect_prodConstraint.tif','Resolution',300);