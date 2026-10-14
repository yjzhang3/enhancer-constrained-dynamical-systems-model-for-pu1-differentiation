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

%% Compute solved production rate at pu1 range
pu1_totest_conc = 0:0.001:ss_ylim;
basal_sol = eval(S(1).basal);
vmax_sol =  eval(S(1).vmax);
khalf_sol =  eval(S(1).khalf);
N_sol = eval(S(1).N);
Hill_pars = [basal_sol,vmax_sol,N_sol,khalf_sol];
% save("prodConstraint_hillpars.mat",'Hill_pars');

anyn_Hill = @(basal,vmax,khalf,N,x) basal+(vmax-basal)*(x^N/(x^N+khalf^N));
prod_totest_conc = zeros(size(pu1_totest_conc));
for i = 1:length(pu1_totest_conc)
    prod_totest_conc(i) = anyn_Hill(basal_sol,vmax_sol,khalf_sol,N_sol,pu1_totest_conc(i));
end

%% plot cell speciifc
% 
close all
f = figure('Position',[81,64,1286,996]);
axis tight


for ct = 1:4
subplot(2,2,ct)

% data uncertainty in ratio
switch ct
    case 1
        ct_patch = color_patch_pu1_bcell_overall;
        ct_color = patch_bcell_color;
        ct_scale = 'linear';
    case 2
        ct_patch = color_patch_pu1_ProLow_overall;
        ct_color = patch_prolow_color;
        ct_scale = 'linear';
    case 3
        ct_patch = color_patch_pu1_ProHigh_overall;
        ct_color = patch_prohigh_color;
        ct_scale = 'linear';
    otherwise
        ct_patch = zeros(4,2);
        ct_scale = 'linear';
end

patch('Faces',f1,'Vertices',ct_patch,'FaceColor',ct_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on

% removal rate and uncerainty
if ct == 4
% plot(pu1_totest_conc,1*pu1_totest_conc,'LineWidth',lw*0.8,'Color',remv_color,'Marker',remv_celltype_marker{4},'MarkerEdgeColor',remv_color,'MarkerFaceColor','none','MarkerSize',ms)
plot(pu1_totest_conc,1*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{4})
hold on
elseif ct == 3
    plot_range_growthrate_patch(pu1_totest_conc,remRate_constant_ratio_all(ct,3)*pu1_totest_conc,remRate_constant_ratio_all(ct,1)*pu1_totest_conc,hex2rgb(remv_color{ct}))
    hold on
    plot(pu1_totest_conc,remRate_constant_ratio_all(ct,1)*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{ct})
    hold on
else
    plot_range_growthrate_patch(pu1_totest_conc,remRate_constant_ratio_all(ct,3)*pu1_totest_conc,remRate_constant_ratio_all(ct,1)*pu1_totest_conc,hex2rgb(remv_color{ct}))
    hold on
    plot(pu1_totest_conc,remRate_constant_ratio_all(ct,2)*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{ct})
    hold on
end

% overlay monotonic production rate found by sovling four equations
plot(pu1_totest_conc,prod_totest_conc,'LineWidth',lw,'Color',prod_color)

% final touch
xlim([0 ss_ylim])
ylim([0 prod_ylim])
xlabel('\it{[pu]}')
ylabel('Norm. rate')

if ct < 4
legend('','','Remv. rate','Prod. rate','Location','northwest','fontsize',fs_small)
else
legend('','Remv. rate','Prod. rate','Location','northwest','fontsize',fs_small)
end

ax = gca;
ax.Box = 'on';
ax.FontSize = fs;
ax.XScale = ct_scale;
ax.YScale = ct_scale;

switch ct_scale
    case 'log'
    xlim([10^(-3) 1])
    ylim([10^(-3) 1])
    otherwise
    xlim([0 1.2])
    ylim([0 2.5])
end

end

exportgraphics(f,'PU1_alone_Hill_combined_fail.tif','Resolution',300);

%% plot two cell speciifc for main figure
% 
close all
%f = figure('Position',[81,64,1286,996]);
axis tight


for ct = [1,3]
f = figure('Position',[100 100 800 600]);

% data uncertainty in ratio
switch ct
    case 1
        ct_patch = color_patch_pu1_bcell_overall;
        ct_color = patch_bcell_color;
        ct_scale = 'linear';
    case 2
        ct_patch = color_patch_pu1_ProLow_overall;
        ct_color = patch_prolow_color;
        ct_scale = 'linear';
    case 3
        ct_patch = color_patch_pu1_ProHigh_overall;
        ct_color = patch_prohigh_color;
        ct_scale = 'linear';
    otherwise
        ct_patch = zeros(4,2);
        ct_scale = 'linear';
end

patch('Faces',f1,'Vertices',ct_patch,'FaceColor',ct_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on

% removal rate and uncerainty
if ct == 4
% plot(pu1_totest_conc,1*pu1_totest_conc,'LineWidth',lw*0.8,'Color',remv_color,'Marker',remv_celltype_marker{4},'MarkerEdgeColor',remv_color,'MarkerFaceColor','none','MarkerSize',ms)
plot(pu1_totest_conc,1*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{4})
hold on
elseif ct == 3
    plot_range_growthrate_patch(pu1_totest_conc,remRate_constant_ratio_all(ct,3)*pu1_totest_conc,remRate_constant_ratio_all(ct,1)*pu1_totest_conc,hex2rgb(remv_color{ct}))
    hold on
    plot(pu1_totest_conc,remRate_constant_ratio_all(ct,1)*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{ct})
    hold on
else
    plot_range_growthrate_patch(pu1_totest_conc,remRate_constant_ratio_all(ct,3)*pu1_totest_conc,remRate_constant_ratio_all(ct,1)*pu1_totest_conc,hex2rgb(remv_color{ct}))
    hold on
    plot(pu1_totest_conc,remRate_constant_ratio_all(ct,2)*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{ct})
    hold on
end

% overlay monotonic production rate found by sovling four equations
plot(pu1_totest_conc,prod_totest_conc,'LineWidth',lw,'Color',prod_color)

% final touch
xlim([0 ss_ylim])
ylim([0 prod_ylim])
xlabel('\it{[pu]}')
ylabel('Norm. rate')

if ct < 4
legend('','','Remv. rate','Prod. rate','Location','northwest','fontsize',fs_small)
else
legend('','Remv. rate','Prod. rate','Location','northwest','fontsize',fs_small)
end

ax = gca;
ax.Box = 'on';
ax.FontSize = fs;
ax.XScale = ct_scale;
ax.YScale = ct_scale;

switch ct_scale
    case 'log'
    xlim([10^(-3) 1])
    ylim([10^(-3) 1])
    otherwise
    xlim([0 1.2])
    ylim([0 2.5])
end

exportgraphics(f,sprintf('PU1_alone_Hill_combined_fail_%.0f.tif',ct),'Resolution',300);
end



%% find intersections of prod and each remv curve
clc
pred_pu1_SS = zeros(3,3); % there are three intersections for each ct-specific remv rate. Each row is a cell ypte; each column is an intersection
pred_pu1_prod = zeros(3,3);

for ct = 1:3
    % if ct == 4
    %     ind = find(abs(1*pu1_totest_conc-prod_totest_conc)<0.001);
    %     concise_ind = extract_intersection_index(ind,40);
    %     pred_pu1_SS(:,ct) = concise_ind;
    %     pred_pu1_prod(:,ct) = prod_totest_conc(concise_ind);
    % else
    ind = find(abs(remRate_constant_ratio_all(ct,2)*pu1_totest_conc-prod_totest_conc)<0.001);
    concise_ind = extract_intersection_index(ind,40);
    pred_pu1_SS(ct,:) = pu1_totest_conc(concise_ind);
    pred_pu1_prod(ct,:) = prod_totest_conc(concise_ind);
    % end
end

%% present results
% show two stable states for each ct
% show_results_real_vs_pred_v3_log_twoSSperCell(pred_pu1_SS(:,[1,3]),...
%     pred_pu1_prod(:,[1,3]),'PU1_alone_noData')

% show only one stable for each ct
show_results_real_vs_pred_v3(pred_pu1_SS(:,[1]),...
    pred_pu1_prod(:,[1]),'PU1_alone_noData','log')