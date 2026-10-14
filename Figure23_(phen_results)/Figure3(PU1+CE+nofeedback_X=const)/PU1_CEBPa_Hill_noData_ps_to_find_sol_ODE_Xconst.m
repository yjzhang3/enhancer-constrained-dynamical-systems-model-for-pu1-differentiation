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
fs = 25;
fs_small = 18;
patch_bcell_color = '#FF5800';
patch_prolow_color = '#BFBFBF';
patch_prohigh_color = 'k';

tau_Bcell = [17.7,25.5,41.2];
tau_ProLow = [19.9,28.8,47.8];
tau_ProHigh = [9.47,11.6,15.8];
D_pu1 = 1/70;
tau_Mac = 40;

%% optimize
clc
nvars = 8;
% fun = @obj_twoHill;
fun = @obj_oneHill_constX_ODE;
lb = zeros(1,nvars)+0;
ub = zeros(1,nvars)+1;

lb([1]) = 0.0001; % basal
ub([1]) = 0.2;

lb([2]) = 0.5; % vmax
ub([2]) = 1;

lb([3]) = 1; % Hill coefficient
ub([3]) = 10;

lb([5:8]) = 0;
ub([5:8]) = 1; % cell type-specific constant
fn = 'ps_ODE_Xconstant_oneHill_Xconst.txt';

% for ff = 1:1
%    [x,fval]  = particleswarm(fun,nvars,lb,ub,optimoptions("particleswarm",'Display','iter'));
%    create_parameter_file(fn,x,fval,'a+')
% end

%% read best data and check if it's true stable state
close all
clc
% fn = 'ps_to_find_f12gHill_sol_Slope_wFeedback.txt';
A = readmatrix(fn);
[errmin,I] = min(A(:,end));
% I = length(A(:,end));
best_pars = A(I,1:end-1);
best_p = best_pars;

% do this outside of objective so it's faster
p_Hill = best_p(1:4);
g_rate_Bcell = best_p(5);
g_rate_ProLow = best_p(6);
g_rate_ProHigh = best_p(7);
g_rate_mac = best_p(8);

% run objective function with best pars
[err,tB,tPL,tPH,tmac,...
    Bcell_over_time,PL_over_time,PH_over_time,Mac_over_time,...
    pred_pu1_ss,pred_pu1_prod] = obj_oneHill_constX_ODE(best_p);


figure();
plot(tB,Bcell_over_time,'LineWidth',3,'Color',remv_color{1})
hold on
plot(tPL,PL_over_time,'LineWidth',3,'Color',remv_color{2})
hold on
plot(tPH,PH_over_time,'LineWidth',3,'Color',remv_color{3})
hold on
plot(tmac,Mac_over_time,'LineWidth',3,'Color',remv_color{4})

xlabel('time')
ylabel('\it{[pu]}')

set(gca,'FontSize',fs)

%% plot prod vs. remv rate to double check (I'm not sure if that's viable anymore because now
% I not only need to give some range of pu1, but also some range of TF X?
% but I can reduce the system to 1D now that my optimization is done
clc

basal = p_Hill(1);
vmax = p_Hill(2) ;
N = p_Hill(3);
khalf = p_Hill(4);
% g_basal = p_Hill(5);
% g_vmax = p_Hill(6);
% g_N = p_Hill(7);
% g_khalf = p_Hill(8);

% a range of Pu1
pu1_totest_conc = 0:0.0001:1.2;
tau_cellType = [tau_Bcell(2),tau_ProLow(2),tau_ProHigh(3),tau_Mac];
pred_g_ct = [g_rate_Bcell,g_rate_ProLow,g_rate_ProHigh,g_rate_mac];

pu_prod_totest_range = zeros(length(tau_cellType),length(pu1_totest_conc));

production_rate = @(x,g) (basal+(vmax-basal)*(x^N/(x^N+khalf^N)))*g;
% production_rate_gonly = @(w) g_basal+(g_vmax-g_basal)*(w^g_N/(w^g_N+g_khalf^g_N));

for ct = 1:4
    pred_g = pred_g_ct(ct);
    for pp = 1:length(pu1_totest_conc)
        x = pu1_totest_conc(pp);
        pu_prod_totest_range(ct,pp) = production_rate(x,pred_g);
    end
end


%% plot cell specific
% Range of [PU1] concentration in which we find intersection (color patch)
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
plot(pu1_totest_conc,pu_prod_totest_range(ct,:),'LineWidth',lw,'Color',prod_color)

% axis label
xlabel('\it{[pu]}')
% ylabel('Normalized rate')

% legend
if ct < 4
legend('','Removal rate','','Production rate','Location','northwest','FontSize',fs_small)
else
legend('','Removal rate','Production rate','Location','northwest','FontSize',fs_small)
end

% plot(pu1_totest_conc(ss_ind(ct)),prod_totest_range(ct,ss_ind(ct)),'LineStyle','none','Marker','o')
    
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
    ylim([10^(-2) 1])
    otherwise
    xlim([0 1.2])
    ylim([0 2.5])  
end

end
exportgraphics(f,'PU1_nofeedback_Xconst.tif','Resolution',300);

%% plot the mysterious TF (X) concentration as a function of PU1, and production rate of X as a function of PU1
f = figure('Position',[100 100 800 600]);
axis tight


% colororder({'b','m'})

% yyaxis left
plot([conc_ratio_all(:,2);1]',pred_g_ct(1:4),'LineWidth',lw,'Color',prod_color)
hold on
ylabel('X')
xlabel('{\it[pu]}');
ylim([0 2])
conc_arr = [conc_ratio_all(:,2);1];
dot_colors = {'#FF5800','#808080','k','#2626B2'};
for i = 1:4
    plot(conc_arr(i),pred_g_ct(i),'LineStyle','none','Marker','o','MarkerFaceColor',dot_colors{i},'MarkerEdgeColor','none','MarkerSize',ms);
    hold on
end

% yyaxis right
% hold on
% plot(conc_ratio_all(:,2)',w_prod_totest_range(1:3),'LineWidth',lw,'Marker','o','MarkerFaceColor','m','MarkerEdgeColor','m')
% ylabel('{\it [x]} prod. rate')
% ylim([0 1])

set(gca,'FontSize',fs)

legend({'','B cell','Pro low','Pro high','Macrophage'}, ...
       'Location', 'NorthWest','FontSize',fs_small); 
exportgraphics(f,'learned_from_nofeedback.tif','Resolution',300);


%% present results in a readable way
show_results_real_vs_pred_v3(pred_pu1_ss,pred_pu1_prod,'PU1_CEBPa_Hill','log')