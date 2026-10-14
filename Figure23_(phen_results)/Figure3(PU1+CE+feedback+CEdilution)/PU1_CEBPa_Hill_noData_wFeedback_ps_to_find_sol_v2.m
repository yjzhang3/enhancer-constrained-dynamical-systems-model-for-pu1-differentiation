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

% tau_Mac = [37.8,40,56.4];
tau_Bcell = [17.7,25.5,41.2];
tau_ProLow = [19.9,28.8,47.8];
tau_ProHigh = [9.47,11.6,15.8];
D_pu1 = 1/70;
tau_Mac = 40;

%% optimize
clc
nvars = 13;
% fun = @obj_twoHill;
fun = @obj_twoHill_wSlope_wFeedback;
lb = zeros(1,nvars)+0.00001;
ub = zeros(1,nvars)+5;
ub([1,5,9]) = 0.5;
lb([3,7,11]) = 1;
ub([3,7,11]) = 10;
ub([2,6,10]) = 1.2;
ub(13) = 0.1; % degradation rate of X restricted to be very small; can fit still be good?
% initial_point = [0.176325,0.920545,9.862062,0.260102,0.320784,1.122661,5.189452,0.393969,0.168374,0.972478,9.735012,0.096231,4.977935];

fn = 'ps_to_find_f12gHill_sol_Slope_wFeedback_negResults.txt';

% for ff = 1:10
%    [x,fval]  = particleswarm(fun,nvars,lb,ub,optimoptions('particleswarm','Display','iter'));
%        % ,'InitialPoints',initial_point));
%    % [x,fval]  = ga(fun,nvars,[],[],[],[],lb,ub,[],optimoptions('ga','Display','iter'));
%    create_parameter_file(fn,x,fval,'a+')
% end

%% read best data and check if it's true stable state
close all
clc
A = readmatrix(fn);
[errmin,I] = min(A(:,end));
% I = length(A(:,end));
best_pars = A(I,1:end-1);
best_p = best_pars;

% do this outside of objective so it's faster
pars_f1_hill = best_p(1:4); % autoregulation
pars_f2_hill = best_p(5:8); % PU1 regulated by X
pars_g_hill = best_p(9:12); % X regulated by PU1
D_X = best_p(13); % natural degradation rate of TF X

% run objective function with best pars
[final_err,pred_pu1_ss,pred_X_ss,pred_pu1_prod,...
    tb,tpl,tph,tm,...
    bcell_over_time,prolow_over_time,prohigh_over_time,mac_over_time] = obj_twoHill_wSlope_wFeedback(best_p);

% see if my event is satisfied
t = 0;
% bcell_event = ODE_func_cellType_myEvents(t,[pred_pu1_ss(1),pred_X_ss(1)],pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,tau_Bcell(2));
% prolow_event = ODE_func_cellType_myEvents(t,[pred_pu1_ss(2),pred_X_ss(2)],pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,tau_ProLow(2));
% prohigh_event = ODE_func_cellType_myEvents(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,tau_ProHigh(2));
% mac_event = ODE_func_cellType_myEvents(t,y,pars_f1_hill,pars_f2_hill,pars_g_hill,D_X,tau_Mac(2));


figure();
plot(tb,bcell_over_time,'LineWidth',3,'Color',remv_color{1})
hold on
plot(tpl,prolow_over_time,'LineWidth',3,'Color',remv_color{2})
hold on
plot(tph,prohigh_over_time,'LineWidth',3,'Color',remv_color{3})
hold on
plot(tm,mac_over_time,'LineWidth',3,'Color',remv_color{4})

xlabel('time')
ylabel('\it{[pu]}')

set(gca,'FontSize',fs)

%% plot prod vs. remv rate to double check (I'm not sure if that's viable anymore because now
% I not only need to give some range of pu1, but also some range of TF X?
% but I can reduce the system to 1D now that my optimization is done
clc

basal = pars_f1_hill(1);
vmax = pars_f1_hill(2) ;
N = pars_f1_hill(3);
khalf = pars_f1_hill(4);

basal_2 = pars_f2_hill(1);
vmax_2 = pars_f2_hill(2) ;
N_2 = pars_f2_hill(3);
khalf_2 = pars_f2_hill(4);

g_basal= pars_g_hill(1);
g_vmax = pars_g_hill(2);
g_N = pars_g_hill(3);
g_khalf = pars_g_hill(4);

% a range of Pu1
pu1_totest_conc = 0:0.00001:1.2;
tau_cellType = [tau_Bcell(2),tau_ProLow(2),tau_ProHigh(3),tau_Mac];
prod_totest_range = zeros(length(tau_cellType),length(pu1_totest_conc));
ss_ind = zeros(length(tau_cellType),1);

for ct = 1:4
    norm_remv_rate_const_PU1 = (log(2)/tau_cellType(ct)+D_pu1)/(log(2)/tau_Mac+D_pu1);
    norm_remv_rate_const_X = (log(2)/tau_cellType(ct)+D_X)/(log(2)/tau_Mac+D_X);
    for pp = 1:length(pu1_totest_conc)
        y = zeros(1,2);
        y(1) = pu1_totest_conc(pp);
        y(2) = (g_basal+(g_vmax-g_basal)*(y(1)^g_N/(y(1)^g_N+g_khalf^g_N)))/norm_remv_rate_const_X;

        prod_totest_range(ct,pp) = (basal+(vmax-basal)*(y(1)^N/(y(1)^N+khalf^N)))*...
        (basal_2+(vmax_2-basal_2)*(y(2)^N_2/(y(2)^N_2+khalf_2^N_2))); 

        
        % if abs(pu1_totest_conc(pp) - pred_pu1_ss(ct)) < 0.0001
        %     % disp(y(2))
        %     % disp(prod_totest_range(ct,pp)-norm_remv_rate_const_PU1*y(1))
        % end
    end

    raw_ind = find(abs(pu1_totest_conc - pred_pu1_ss(ct)) < 0.0001);
    ss_ind(ct) = raw_ind(1);

end

%% Range of [PU1] concentration in which we find intersection (color patch)
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

%% plot cell specific
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
plot(pu1_totest_conc,prod_totest_range(ct,:),'LineWidth',lw,'Color',prod_color)

% axis label
xlabel('\it{[pu]}')
ylabel('Normalized rate')

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
exportgraphics(f,'PU1_CEBPa_Hill_can_work.tif','Resolution',300);

%% present results in a readable way
show_results_real_vs_pred_v3(pred_pu1_ss,pred_pu1_prod,'PU1_CEBPa_Hill','log')