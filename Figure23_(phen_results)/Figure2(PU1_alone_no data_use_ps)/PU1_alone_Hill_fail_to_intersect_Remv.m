% four equatinos can either be let prod = rmv
% doesn't work because 1) prohigh rmv rate > prolow rmv rate, so prod
% intersects prohigh first and then prolow; 2) however, pu1 conc at prohigh
% must be larger than conc at prolow (i.e., intersections must happen for
% prolow first)
% not necessarily!! What if it's not not a sigmoid shape but S-shpaed? To
% strengthen the argument, I'll just use the pars found from prodConstraint
% and plot the production rate here to argue that production intersection
% is incorrect.

% eqn1 = Hill_func_Bcell == remRate_constant_ratio_all(1,2)*conc_ratio_all(1,2);
% eqn2 = Hill_func_Bcell == remRate_constant_ratio_all(2,2)*conc_ratio_all(2,2);
% eqn3 = Hill_func_Bcell == remRate_constant_ratio_all(3,2)*conc_ratio_all(3,2);
% eqn4 = Hill_func_Bcell == 1*1;

% this can be shown by directly plotting four remv curves

%%
clear
clc
close all
load(fullfile(fileparts(pwd), 'Figure1A (intro) and S1 (replicate errorbar)', 'error_prop.mat'));
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
prod_color = [0.00 0.75 0.65];

%%
f = figure('Position',[100 100 800 600]);
axis tight

%%%% removal rate
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
hold on

%%%% overlay monotonic production rate found by sovling four equations
load("prodConstraint_hillpars.mat",'Hill_pars');
basal_sol = Hill_pars(1);
vmax_sol =  Hill_pars(2);
N_sol = Hill_pars(3);
khalf_sol =  Hill_pars(4);


pu1_totest_conc = 0:0.01:1.2;
anyn_Hill = @(basal,vmax,khalf,N,x) basal+(vmax-basal)*(x^N/(x^N+khalf^N));
prod_totest_conc = zeros(size(pu1_totest_conc));
for i = 1:length(pu1_totest_conc)
    prod_totest_conc(i) = anyn_Hill(basal_sol,vmax_sol,khalf_sol,N_sol,pu1_totest_conc(i));
end

plot(pu1_totest_conc,prod_totest_conc,'LineWidth',lw,'Color',prod_color)

xlim([0 1.2])
ylim([0 2.5])

xlabel('\it{[pu]}')
ylabel('Normalized rate')

legend('B cell remv. rate','Pro low remv. rate','Pro high remv. rate','Macrophage remv. rate','Phen. prod. rate','Location','northwest')

set(gca,'FontSize',fs)

exportgraphics(f,'PU1_alone_Hill_fail_to_intersect_Remv.tif','Resolution',300);
