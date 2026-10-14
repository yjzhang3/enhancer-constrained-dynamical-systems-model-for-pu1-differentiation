% use global optimization to find solution such that product of two hill
% equations (concentration at each cell type is different beteween f and
% g!) can produce a production rate intereseting ct-specific growth rate
% and lie in the data constraint
clear
clc
close all
prod_color = [0.00 0.75 0.65];
remv_color = {'#FF5800','#BFBFBF','#808080','#2626B2'}; % or different color for different ct
% remv_celltype_marker = {'+','d','*','s'}; % either different marker for different ct
lw = 3;
ms = 13;
fa = 0.6;
fs = 19;
patch_bcell_color = '#FF5800';
patch_prolow_color = '#BFBFBF';
patch_prohigh_color = 'k';
delete ps_to_find_fgHill_sol_n3.txt

%% optimize
nvars = 13;
% fun = @obj_twoHill;
fun = @obj_twoHill_wSlope;
lb = zeros(1,nvars)+0;
ub = zeros(1,nvars)+10;

% no cooperativity yet
lb(9) = 1;
ub(9) = 1;

% Hill coefficient
lb([3]) = 1;
ub([3]) = 20;

% max of cebpa should be 1 in macrophage (but you don't know the mysterious TF is cebpa,
% so you can't say the TF's maximum happens in macrophage). so don't
% constrain macrophage specifically 
% ub(10:13) = 1;

% maximum f(x) is 1 (pu1 production rate maximum happens in macrophage)
% ub(2) = 1;

fn = sprintf('ps_to_find_fgHill_sol_N%.0f.txt',lb(3));

for ff = 1:10
   [x,fval]  = particleswarm(fun,nvars,lb,ub);
   create_parameter_file(fn,x,fval,'a+')
end

%% read best data
clc
A = readmatrix(fn);
[errmin,I] = min(A(:,end));
% I = 2;
best_pars = A(I,1:end-1);
p = best_pars;

% evaluating function handle is too slow so explicitly initating pars
basal = p(1);
vmax = p(2) ;
N = p(3);
khalf = p(4);
g_basal= p(5);
g_vmax=p(6);
g_N = p(7)  ;
g_khalf = p(8);
coop = p(9);
g_input_Bcell = p(10);
g_input_ProLow = p(11);
g_input_ProHigh = p(12);
g_input_mac = p(13);

pu1_totest_conc = 0:0.001:1.2;
prod_totest_range = zeros(4,length(pu1_totest_conc));
cebpa_totest_conc = p(10:13);
for ct = 1:4
    cebpa = cebpa_totest_conc(ct);
    for pp = 1:length(pu1_totest_conc)
        pu1 = pu1_totest_conc(pp);
        prod_totest_range(ct,pp) = (basal+(vmax-basal)*(pu1^N/(pu1^N+khalf^N)))*...
            (g_basal+(g_vmax-g_basal)*(cebpa^g_N/(cebpa^g_N+g_khalf^g_N)));
    end
end


%% Range of [PU1] concentration in which we find intersection (color patch)
close all
load(fullfile(fileparts(pwd), 'Error propagation and Normalization', 'error_prop.mat'));

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

f = figure('Position',[124,67,1552,1029]);
axis tight

% uncertainty in ratio
ms_adjust_factor = 0.5;
subplot(2,2,1)
patch('Faces',f1,'Vertices',color_patch_pu1_bcell_overall,'FaceColor',patch_bcell_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
% plot(conc_ratio_all(1,2),prodRate_ratio_all(1,2),'LineStyle','none','Marker','s','MarkerFaceColor',patch_bcell_color,'MarkerEdgeColor',patch_bcell_color,'MarkerSize',ms*ms_adjust_factor,'LineWidth',6);

subplot(2,2,2)
patch('Faces',f1,'Vertices',color_patch_pu1_ProLow_overall,'FaceColor',patch_prolow_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
% plot(conc_ratio_all(2,2),prodRate_ratio_all(2,2),'LineStyle','none','Marker','s','MarkerFaceColor',patch_prolow_color,'MarkerEdgeColor',patch_prolow_color,'MarkerSize',ms*ms_adjust_factor,'LineWidth',6);

subplot(2,2,3)
patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',patch_prohigh_color,'EdgeColor','none','FaceAlpha',fa,'LineWidth',lw)
hold on
% plot(conc_ratio_all(3,2),prodRate_ratio_all(3,2),'LineStyle','none','Marker','s','MarkerFaceColor',patch_prohigh_color,'MarkerEdgeColor',patch_prohigh_color,'MarkerSize',ms*ms_adjust_factor,'LineWidth',6);

subplot(2,2,4)
plot(1,1,'LineStyle','none','Marker','pentagram','MarkerEdgeColor',remv_color{4},'MarkerFaceColor',remv_color{4},'MarkerSize',ms)
hold on

cell_name = {'B cell','pro low','pro high'};
% remv and prod rate
for ct = 1:4
    % same color, different mk
    % plot(pu1_totest_conc,remRate_constant_ratio_all(ct,2)*pu1_totest_conc,'LineWidth',lw,'Color',remv_color,'Marker',remv_celltype_marker{ct},'MarkerEdgeColor',remv_color,'MarkerFaceColor','none','MarkerSize',ms)
    subplot(2,2,ct)
    
    % same mk, different color
    if ct == 4
        plot(pu1_totest_conc,1*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{ct})
        hold on
        plot(pu1_totest_conc,prod_totest_range(ct,:),'LineWidth',lw,'Color',prod_color)
        xlabel('\it{[pu]}')
        ylabel('Normalized rate')
        legend('','Remv. rate','Phen. prod. rate','Location','best')
    else
        plot(pu1_totest_conc,remRate_constant_ratio_all(ct,2)*pu1_totest_conc,'LineWidth',lw,'Color',remv_color{ct})
        hold on
        plot(pu1_totest_conc,prod_totest_range(ct,:),'LineWidth',lw,'Color',prod_color)
        xlabel('\it{[pu]}')
        ylabel('Normalized rate')
        legend(strcat('Data-constrained'," ",cell_name{ct}),'Remv. rate','Phen. prod. rate','Location','best')
    end
    
end

ax = findall(gcf, 'Type', 'axes');
set(ax,'Box','on')
xlim(ax,[10^(-3) 1])
ylim(ax(2),[10^(-3) 10^0.5])
ylim(ax([1,3,4]),[10^(-3) 10^(0.3)])
set(ax, 'FontSize', fs)
set(ax,'YScale','log','XScale','log')
exportgraphics(f,'PU1_CEBPa_Hill_can_work.tif','Resolution',300);

%% present results in a readable way
pred_pu1_SS = zeros(1,3);
pred_pu1_prod = zeros(1,3);
for ct = 1:3
    ind = find(abs(prod_totest_range(ct,:)-remRate_constant_ratio_all(ct,2)*pu1_totest_conc)<0.0001);
    % concise_ind = extract_intersection_index(ind,40);
    concise_ind = ind(1);
    pred_pu1_SS(ct) = pu1_totest_conc(concise_ind);
    pred_pu1_prod(ct) = prod_totest_range(ct,concise_ind);
end


show_results_real_vs_pred_v3(pred_pu1_SS,pred_pu1_prod,'PU1_CEBPa_Hill')