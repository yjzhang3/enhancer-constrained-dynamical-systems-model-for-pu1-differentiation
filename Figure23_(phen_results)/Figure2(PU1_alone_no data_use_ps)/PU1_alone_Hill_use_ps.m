%% housekeeping
clear
clc
close all
load('error_prop_0216.mat');
prod_color = [0.00 0.75 0.65];
% remv_color = [0.70 0.40 0.85];
remv_color = {'#FF5800','#BFBFBF','#808080','#2626B2'}; % or different color for different ct
% remv_celltype_marker = {'+','d','*','s'}; % either different marker for different ct
lw = 3;
ms = 15;
fa = 0.6;
fs = 25;
fs_small = 18;
patch_bcell_color = '#FF5800';
patch_prolow_color = '#BFBFBF';
patch_prohigh_color = 'k';
ss_ylim = 1.2;
prod_ylim = 2.5;

%% fit
fun = @(p) obj_ps_find_Hill_pars(p);
nvars = 8;
lb = zeros(1,nvars);ub = zeros(1,nvars);

lb(1:4) = [0.0001,0.2,1,0.0001]; ub(1:4) = [0.2,1,10,2];
lb(5:8) = [conc_ratio_all(1,1),conc_ratio_all(2,1),conc_ratio_all(3,1),1]; 
ub(5:8) = [conc_ratio_all(1,3),conc_ratio_all(2,3),conc_ratio_all(3,3),1];

fn = 'ps_to_find_Hill.txt';

% for ff = 1:1
%    [x,fval]  = particleswarm(fun,nvars,lb,ub,optimoptions('particleswarm','Display','iter'));
%    create_parameter_file(fn,x,fval,'a+')
% end

%% predict
close all
clc
A = readmatrix(fn);
[errmin,I] = min(A(:,end));
% I = length(A(:,end));
best_pars = A(I,1:end-1);

hill_pars = best_pars(1:4);
ct_pu = best_pars(5:8);

% a range of Pu1
pu1_totest_conc = 0:0.001:1.2;
prod_totest_range = zeros(1,length(pu1_totest_conc));
Hill_func = @(hill_pars,x) hill_pars(1)+(hill_pars(2)-hill_pars(1))*(x^hill_pars(3)/(x^hill_pars(3)+hill_pars(4)^hill_pars(3)));

for puu = 1:length(pu1_totest_conc)
    prod_totest_range(puu) = Hill_func(hill_pars,pu1_totest_conc(puu));
end

%% plot cell specific

close all
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


close all
f = figure('Position',[81,64,1286,996]);
axis tight

for ct = 1:4
subplot(2,2,ct)

% uncertainty in ratio
switch ct
    case 1
        ct_patch = color_patch_pu1_bcell_overall;
        ct_color = patch_bcell_color;
        ct_remv = remRate_constant_ratio_all(ct,2)*pu1_totest_conc;
        ct_pltoptions = 'linear';
    case 2
        ct_patch = color_patch_pu1_ProLow_overall;
        ct_color = patch_prolow_color;
        ct_remv = remRate_constant_ratio_all(ct,2)*pu1_totest_conc;
        ct_pltoptions = 'linear';
    case 3
        ct_patch = color_patch_pu1_ProHigh_overall;
        ct_color = patch_prohigh_color;
        ct_remv = remRate_constant_ratio_all(ct,1)*pu1_totest_conc;
        ct_pltoptions = 'linear';
    otherwise
        ct_patch = zeros(4,2);
        ct_remv =  1*pu1_totest_conc;
        ct_pltoptions = 'linear';
end

patch('Faces',f1,'Vertices',ct_patch,'FaceColor',ct_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on

% removal rate
plot(pu1_totest_conc,ct_remv,'LineWidth',lw,'Color',remv_color{ct})
hold on

% removal rate uncertainty
if ct < 4
max_removal_rate = remRate_constant_ratio_all(ct,3)*pu1_totest_conc;
min_removal_rate = remRate_constant_ratio_all(ct,1)*pu1_totest_conc;
plot_range_growthrate_patch(pu1_totest_conc,max_removal_rate,min_removal_rate,hex2rgb(remv_color{ct}))
hold on
end

% production 
plot(pu1_totest_conc,prod_totest_range,'LineWidth',lw,'Color',prod_color)

% axis label
% xlabel('\it{[pu]}')
% ylabel('Synthesis rate')

% legend
if ct < 4
legend('','Removal rate','','Production rate','Location','northwest','FontSize',fs_small);
else
legend('','Removal rate','Production rate','Location','northwest','FontSize',fs_small);
end
% final touch
xlabel('\it{[pu]}')
ylabel('Rate')

ax = gca;
ax.Box = 'on';
ax.FontSize = fs;
ax.XScale = ct_pltoptions;
ax.YScale = ct_pltoptions;

switch ct_pltoptions
    case 'log'
    xlim([0 1.2])
    ylim([10^(-4) 10])
    otherwise
    xlim([0 1.2])
    ylim([0 2.5])  
end

end
exportgraphics(f,'PU1_alone_Hill_cant_work_use_ps.tif','Resolution',300);

%% plot two cell speciifc for main figure 


close all
% f = figure('Position',[81,64,1286,996]);


for ct = 1:4
f = figure('Position',[100 100 800 600]);
axis tight

% uncertainty in ratio
switch ct
    case 1
        ct_patch = color_patch_pu1_bcell_overall;
        ct_color = patch_bcell_color;
        ct_remv = remRate_constant_ratio_all(ct,2)*pu1_totest_conc;
        ct_pltoptions = 'linear';
    case 2
        ct_patch = color_patch_pu1_ProLow_overall;
        ct_color = patch_prolow_color;
        ct_remv = remRate_constant_ratio_all(ct,2)*pu1_totest_conc;
        ct_pltoptions = 'linear';
    case 3
        ct_patch = color_patch_pu1_ProHigh_overall;
        ct_color = patch_prohigh_color;
        ct_remv = remRate_constant_ratio_all(ct,1)*pu1_totest_conc;
        ct_pltoptions = 'linear';
    otherwise
        ct_patch = zeros(4,2);
        ct_remv =  1*pu1_totest_conc;
        ct_pltoptions = 'linear';
end

patch('Faces',f1,'Vertices',ct_patch,'FaceColor',ct_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on

% removal rate
plot(pu1_totest_conc,ct_remv,'LineWidth',lw,'Color',remv_color{ct})
hold on

% removal rate uncertainty
if ct < 4
max_removal_rate = remRate_constant_ratio_all(ct,3)*pu1_totest_conc;
min_removal_rate = remRate_constant_ratio_all(ct,1)*pu1_totest_conc;
plot_range_growthrate_patch(pu1_totest_conc,max_removal_rate,min_removal_rate,hex2rgb(remv_color{ct}))
hold on
end

% production 
plot(pu1_totest_conc,prod_totest_range,'LineWidth',lw,'Color',prod_color)

% % axis label
% xlabel('\it{[pu]}')
% ylabel('Synthesis rate')

% legend
if ct < 4
legend('','Removal rate','','Production rate','Location','northwest','FontSize',fs_small);
else
legend('','Removal rate','Production rate','Location','northwest','FontSize',fs_small);
end
% final touch
xlabel('\it{[pu]}')
ylabel('Rate')

ax = gca;
ax.Box = 'on';
ax.FontSize = fs;
ax.XScale = ct_pltoptions;
ax.YScale = ct_pltoptions;

switch ct_pltoptions
    case 'log'
    xlim([0 1.2])
    ylim([10^(-4) 10])
    otherwise
    xlim([0 1.2])
    ylim([0 2.5])  
end
exportgraphics(f,sprintf('PU1_alone_Hill_cant_work_use_ps_%.0f.tif',ct),'Resolution',300);
end

%% find intersections of prod and each remv curve
clc
pred_pu1_SS = zeros(3,3); % there are three intersections for each ct-specific remv rate. 
% Each row is a cell ypte; each column is an intersection
pred_pu1_prod = zeros(3,3);

for ct = 1:3
    if ct == 3
        ind = find(abs(remRate_constant_ratio_all(ct,1)*pu1_totest_conc-prod_totest_range)<0.001);
        concise_ind = extract_intersection_index(ind,40);
        pred_pu1_SS(ct,:) = pu1_totest_conc(concise_ind);
        pred_pu1_prod(ct,:) = prod_totest_range(concise_ind);
    else
    ind = find(abs(remRate_constant_ratio_all(ct,2)*pu1_totest_conc-prod_totest_range)<0.001);
    concise_ind = extract_intersection_index(ind,40);
    pred_pu1_SS(ct,:) = pu1_totest_conc(concise_ind);
    pred_pu1_prod(ct,:) = prod_totest_range(concise_ind);
    end
end

%% present results
% show two stable states for each ct
% show_results_real_vs_pred_v3_log_twoSSperCell(pred_pu1_SS(:,[1,3]),...
%     pred_pu1_prod(:,[1,3]),'PU1_alone_noData')

% show only one stable for each ct
show_results_real_vs_pred_v3(pred_pu1_SS(:,1),...
    pred_pu1_prod(:,1),'PU1_alone_noData','log')
