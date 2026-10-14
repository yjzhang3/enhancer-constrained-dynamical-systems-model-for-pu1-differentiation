function ODE_prediction_cebpa_PU1CEfeedback_1D(objFuncHandle,selec,best_pars,nuc_destab_location)
%% housekeeping
load(fullfile(parentdir(pwd, 1), 'Error propagation and Normalization','error_prop.mat'));

% subdir_fig = strcat('ODEoutputs_figs_',erase(final_fn,'.txt')); % too
% long!!!
subdir_fig = 'ODEoutputs_figs'; 
% subdir_mat = 'ODEoutputs_mats';

% Ensure subdirectory exists (creates if not)
if ~exist(subdir_fig, 'dir')
    mkdir(subdir_fig);
end

% if ~exist(subdir_mat, 'dir')
%     mkdir(subdir_mat);
% end

RNAPFuncHandle = @RNAP_occ_prob_cebpa_general;
%% plotting housekeeping
linestyle_name = {'-','--',':','-o','->','-v','-square'};
lw = 3; % linewidth
ms = 12; % marker size
fs = 25; % fontsize
lg_fs = 16;
fa = 0.6; % facealpha of constraint patch (not growth rate patch)
% color for constraint patch
bcellcolor = '#FF5800';
prolowcolor = '#BFBFBF';
prohighcolor = 'k';
maccolor = '#2626B2';

% color and style for removal rate linear line
% removalcolor = [0.70 0.40 0.85];
removalcolor = {'#FF5800','#BFBFBF','#808080','#2626B2'};
prodcolor = [0.00 0.75 0.65];
removalls = '-';

% color of intersection
color_mk_intersect = '#E03B3B';

% upper limit of plot
xlim_conc = 1.2;
ylim_prod_rate = 2.5;

%% load best pars
% if ~isfile(final_fn)
%     disp('file does not exist')
%     good_pars = 0;
%     D7_prob_bcell = 0;
%     minerr = NaN;
%     return
% else
%     n_datap = length(selec); % number of data points to fit
%     A = readmatrix(final_fn);
% end
% 
% arri = A(:,end);
% sorted_arri = sort(arri);
% minerr = sorted_arri(1); % pick the smallest or second smallest?
% minind = find(arri==minerr);
% best_pars = A(minind,1:end-1);

% nvars = length(A(1,:))-1;
n_datap = length(selec);
[final_err,pu1_ss_conc,cebpa_ss_conc,pu1_ss_prod,...
    real_d_final,sim_d_final,D7_valid_sublist,energy_table,...
    Bcell_over_time,ProLow_over_time,ProHigh_over_time,tB,tPL,tPH] = ...
    objFuncHandle(selec,best_pars,nuc_destab_location,conc_ratio_all,prodRate_ratio_all);
disp(final_err)

%% plot cebpa concentration
figure();
bar({'B cell','Pro low','Pro high'},cebpa_ss_conc)


%% plot prediction vs. experimental gene expression

load('Leddin_real_d_std.mat','real_d_std')
% real_d_std = [real_d_std,10];
all_data_label = {'D2/D1','D3/D1','D4/D1','D5/D1','D6/D1','D7/D1','F2/E1','F3/E1','F4/E1','E2/E1','E3/E1'};

f=figure();
plot(1:n_datap,sim_d_final,'LineStyle','none','Marker','d','MarkerSize',ms,'Color',prodcolor)
hold on
errorbar(1:n_datap,real_d_final,real_d_std(selec)./2556.3,'LineStyle','none','Marker','o','MarkerSize',ms,'Color','k')
set(gca, 'YScale', 'linear')

xticks([1:n_datap])
xticklabels(all_data_label(selec))
ylabel('Ratio in gene expression')

set(gca,'FontSize',fs)
destab_id = num2str(nuc_destab_location);
final_id = strrep(destab_id," ","");

% legend(strcat('Model','(',final_id,')'),'Experiment','Location','Northwest')
legend('Model','Experiment','Location','Northwest','FontSize',lg_fs)

% title(strrep(orig_name, '_', ' '))
% title(final_id)
% notxt = erase(final_fn,".txt"); 
saveas(f,fullfile(subdir_fig,'Gprediction.tif'))

%% plot dynamics (pu vs. time) of cell-type specific model
f = figure('Position',[100 100 800 600]);
axis tight
plot(tB,Bcell_over_time(:,1),'LineWidth',lw,'Color',bcellcolor)
hold on
% plot(tB,Bcell_over_time(:,2),'LineWidth',lw,'Color',bcellcolor,'LineStyle','--')
% hold on
plot(tPL,ProLow_over_time(:,1),'LineWidth',lw,'Color',prolowcolor)
hold on
% plot(tPL,ProLow_over_time(:,2),'LineWidth',lw,'Color',prolowcolor,'LineStyle','--')
% hold on
plot(tPH,ProHigh_over_time(:,1),'LineWidth',lw,'Color',prohighcolor)
hold on
% plot(tPH,ProHigh_over_time(:,2),'LineWidth',lw,'Color',prohighcolor,'LineStyle','--')
% 
xlabel('t')
ylabel('\it{[pu]}')
ylim([0 1])
legend('B cell','','pro low','','pro high','','FontSize',lg_fs)
set(gca,'FontSize',fs)
exportgraphics(f,fullfile(subdir_fig,'TimeDynamics.tif'),'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot rate vs. pu by prod and degradation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Range of [PU1] concentration and synthesis rate in which we find intersection (color patch)
pu1_bcell_conc_range = conc_ratio_all(1,[1,3]); % skip the median of ratio because I'm plotting patch
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

%% range of growth rate
rem_rate_constant_Bcell = remRate_constant_ratio_all(1,2);
rem_rate_constant_Bcell_Q1 = remRate_constant_ratio_all(1,1);
rem_rate_constant_Bcell_Q3 = remRate_constant_ratio_all(1,3);

rem_rate_constant_ProLow = remRate_constant_ratio_all(2,2);
rem_rate_constant_ProLow_Q1 = remRate_constant_ratio_all(2,1);
rem_rate_constant_ProLow_Q3 = remRate_constant_ratio_all(2,3);

rem_rate_constant_ProHigh = remRate_constant_ratio_all(3,2);
rem_rate_constant_ProHigh_Q1 = remRate_constant_ratio_all(3,1);
rem_rate_constant_ProHigh_Q3 = remRate_constant_ratio_all(3,3);

rem_rate_constant_Mac = 1;

%% all possible pu1 level
pu1_totest_range = [0:0.0001:xlim_conc];
D7_ss_prod_rate = sim_d_final(6);
f = figure('Position',[81,64,1286,996]);
axis tight

%% bcell
subplot(2,2,1)
plot_option = 'log';

% %%%%%%% plot overaall color patch
patch('Faces',f1,'Vertices',color_patch_pu1_bcell_overall,'FaceColor',bcellcolor,'EdgeColor','none','FaceAlpha',fa*0.75)
hold on
lgd = "Data-constrained B cell";
hold on

% %%%%%%% generate growth curve
% tau_Bcell = 25.5;
% tau_Bcell_Q1 = 17.7;
% tau_Bcell_Q3 = 41.2;
min_removal_rate_slope = rem_rate_constant_Bcell_Q1;
max_removal_rate_slope = rem_rate_constant_Bcell_Q3;

removal_range = rem_rate_constant_Bcell*pu1_totest_range;
min_removal_rate = min_removal_rate_slope*pu1_totest_range;
max_removal_rate = max_removal_rate_slope*pu1_totest_range;

plot(pu1_totest_range,removal_range,'LineWidth',lw,'Color',removalcolor{1},'LineStyle',removalls)
hold on
plot_range_growthrate_patch(pu1_totest_range,max_removal_rate,min_removal_rate,hex2rgb(removalcolor{1}))
hold on

% %%%%%%% plot production rate
cebp_cell_type_specific_conc_bcell = [cebpa_ss_conc(1)];
bcell_all_possbility = cebp_cell_type_specific_conc_bcell;
bcell_all_possibility_legend = cell(length(bcell_all_possbility(:,1)),1);
bcell_all_possibility_linestyle = cell(length(bcell_all_possbility(:,1)),1);
for i = 1:length(bcell_all_possbility(:,1))
    bcell_all_possibility_legend{i} = sprintf('B cell, \\it[ce]\\rm = %.3f',bcell_all_possbility(i,1));
    bcell_all_possibility_linestyle{i} = linestyle_name{i};
end

TR_range = zeros(size(pu1_totest_range));

for m = 1:length(bcell_all_possbility(:,1))
    cebpa = bcell_all_possbility(m,1);
    % runx1 = bcell_all_possbility(m,2);
    for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
        
        D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        
        D7_curr_prod_rate = D7_prob;

        TR_range(i) = D7_curr_prod_rate/D7_ss_prod_rate;
    end
    
    plot(pu1_totest_range,TR_range,bcell_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
    hold on

    cell_name = bcell_all_possibility_legend{m};
    % lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
end

% %%%%%% find actual steady state
% [~,ind] = min(abs(TR_range-removal_range));
% concise_ind = extract_intersection_index(ind,40);
% true_Bcell_pu1_ss = pu1_totest_range(concise_ind);
% true_Bcell_pu1_prod = TR_range(concise_ind);

% %%%%%% mark itnersection
% xline(model_predicted_Bcell_pu1,'Color','r','LineWidth',2)
% yline(D7_prob_bcell/D7_ss_prod_rate,'Color','r','LineWidth',2)
% plot(pu1_ss_conc(1),pu1_ss_prod(1),'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
% text(model_predicted_Bcell_pu1,norm_D7_prob_bcell,"WT PU.1 production in B cell",'FontSize',12,'Color','r','HorizontalAlignment','left')

% %%%%%% together
% lgd = [lgd,"Remv. rate in B cell","",""];
lgd = ["","Remv. rate","","Prod. rate"];
legend(lgd,'location','northwest','FontSize',lg_fs)
xlabel('\it{[pu]}')
ylabel('Normalized Rate')
set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)

if contains(plot_option,'log')
    plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
    plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
else
    plot_ub_x = [0,pu1_totest_range(end)];
    plot_ub_y = [0,ylim_prod_rate];
end
xlim(plot_ub_x)
ylim(plot_ub_y)

% title(strrep( '_', ' '))
% title(final_id)
% saveas(f,fullfile(subdir_fig,strcat('_Bcellonly','_',plot_option,'.tif')))
ax = gca;
ax.Box = 'on';
exportgraphics(f,fullfile(subdir_fig,strcat('_Bcellonly','_',plot_option,'.tif')),'Resolution',300);


%% plot ProLow 
clear lgd
subplot(2,2,2)
% f = figure('Position',[100 100 800 600]);
% axis tight
plot_option = 'log';

% %%%%%%% plot overall patch
patch('Faces',f1,'Vertices',color_patch_pu1_ProLow_overall,'FaceColor',prolowcolor,'EdgeColor',prolowcolor,'FaceAlpha',fa)
hold on
lgd = "Data-constrained pro low";
hold on

% %%%%%%% generate growth curve
% tau_ProLow = 28.8;
% tau_ProLow_Q1 = 19.9;
% tau_ProLow_Q3 = 47.8;

removal_range_ProLow = rem_rate_constant_ProLow*pu1_totest_range;
min_removal_range_ProLow = rem_rate_constant_ProLow_Q1*pu1_totest_range;
max_removal_range_ProLow = rem_rate_constant_ProLow_Q3*pu1_totest_range;


plot(pu1_totest_range,removal_range_ProLow,'LineWidth',lw,'Color',removalcolor{2},'LineStyle',removalls)
hold on
plot_range_growthrate_patch(pu1_totest_range,max_removal_range_ProLow,min_removal_range_ProLow,hex2rgb(removalcolor{2}))
hold on

% %%%%%%% plot produciont rate for LOW
cebp_cell_type_specific_conc_prolow = [cebpa_ss_conc(2)];
prolow_all_possibility = cebp_cell_type_specific_conc_prolow;
prolow_all_possibility_legend = cell(length(prolow_all_possibility(:,1)),1);
prolow_all_possibility_linestyle = cell(length(prolow_all_possibility(:,1)),1);

for i = 1:length(prolow_all_possibility(:,1))
    prolow_all_possibility_legend{i} = sprintf('pro low, \\it[ce]\\rm = %.3f',prolow_all_possibility(i,1));
    prolow_all_possibility_linestyle{i} = linestyle_name{i};
end
TR_range_LOW = zeros(size(pu1_totest_range));

for m = 1:length(prolow_all_possibility(:,1))
    cebpa = prolow_all_possibility(m,1);
    % runx1 = prolow_all_possibility(m,2);
    for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
        
        D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        
        D7_curr_prod_rate = D7_prob;
        TR_range_LOW(i) = D7_curr_prod_rate/D7_ss_prod_rate;

    end
    plot(pu1_totest_range,TR_range_LOW,prolow_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
    hold on

    cell_name = prolow_all_possibility_legend{m};
    lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
end


% %%%%%% find actual steady state
% [~,ind] = mink(abs(TR_range_LOW-removal_range_ProLow),1000);
% % concise_ind = extract_intersection_index(ind,30);
% true_prolow_pu1_ss = pu1_totest_range(ind(end));
% true_prolow_pu1_prod = TR_range_LOW(ind(end));

% %%%%%% mark itnersection
% plot(pu1_ss_conc(2),pu1_ss_prod(2),'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
% text(model_predicted_ProLow_pu1,norm_D7_prob_ProLow,"WT PU.1 production in Pro Low",'FontSize',12,'Color','r','HorizontalAlignment','left')

% %%%%%%% together
% lgd = [lgd,"Remv. rate in pro low","",""];
lgd = ["","Remv. rate","","Prod. rate"];
legend(lgd,'location','northwest','FontSize',lg_fs)

xlabel('\it{[pu]}')
ylabel('Normalized Rate')
set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)

if contains(plot_option,'log')
    plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
    plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
else
    plot_ub_x = [0,pu1_totest_range(end)];
    plot_ub_y = [0,ylim_prod_rate];
end
xlim(plot_ub_x)
ylim(plot_ub_y)

% title(strrep( '_', ' '))
% title(final_id)
% saveas(f,fullfile('outputs_figs',strcat('_ProLow','_',plot_option,'.tif')))
ax = gca;
ax.Box = 'on';
exportgraphics(f,fullfile(subdir_fig,strcat('_ProLow','_',plot_option,'.tif')),'Resolution',300);

%% then ProHigh
clear lgd
subplot(2,2,3)
% f = figure('Position',[100 100 800 600]);
% axis tight
plot_option = 'linear';

% %%%%%%% plot overall patch
patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',prohighcolor,'EdgeColor',prohighcolor,'FaceAlpha',fa)
hold on
lgd = ["Data-constrained pro high"];
hold on

% %%%%%%% generate growth curve
% tau_ProHigh = 11.6;
% tau_ProHigh_Q1 = 9.47;
% tau_ProHigh_Q3 = 15.8;

removal_range_ProHigh = rem_rate_constant_ProHigh*pu1_totest_range;
min_removal_range_ProHigh = rem_rate_constant_ProHigh_Q1*pu1_totest_range;
max_removal_range_ProHigh = rem_rate_constant_ProHigh_Q3*pu1_totest_range;

plot(pu1_totest_range,min_removal_range_ProHigh,'LineWidth',lw,'Color',removalcolor{3},'LineStyle',removalls)
hold on
plot_range_growthrate_patch(pu1_totest_range,max_removal_range_ProHigh,min_removal_range_ProHigh,hex2rgb(removalcolor{3}))


%%%%%%% plot produciont rate for HIGH
cebp_cell_type_specific_conc_prohigh = [cebpa_ss_conc(3)];
prohigh_all_possibility = cebp_cell_type_specific_conc_prohigh;
prohigh_all_possibility_legend = cell(length(prohigh_all_possibility(:,1)),1);
prohigh_all_possibility_linestyle = cell(length(prohigh_all_possibility(:,1)),1);

for i = 1:length(prohigh_all_possibility(:,1))
    prohigh_all_possibility_legend{i} = sprintf('pro high, \\it[ce]\\rm = %.3f',prohigh_all_possibility(i,1));
    prohigh_all_possibility_linestyle{i} = linestyle_name{i};
end

TR_range_HIGH = zeros(size(pu1_totest_range));
for m = 1:length(prohigh_all_possibility(:,1))
    cebpa = prohigh_all_possibility(m,1);
    % runx1 = prohigh_all_possibility(m,2);
    for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
        
        D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        
        D7_curr_prod_rate = D7_prob;
        TR_range_HIGH(i) = D7_curr_prod_rate/D7_ss_prod_rate;

    end
    plot(pu1_totest_range,TR_range_HIGH,prohigh_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
    hold on

    cell_name = prohigh_all_possibility_legend{m};
    lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
end



% %%%%%% find actual steady state
% [~,ind] = min(abs(TR_range_HIGH-removal_range_ProHigh));
% concise_ind = extract_intersection_index(ind,40);
% true_prohigh_pu1_ss = pu1_totest_range(concise_ind);
% true_prohigh_pu1_prod = TR_range_HIGH(concise_ind);

% %%%%%% mark itnersection
% plot(pu1_ss_conc(3),pu1_ss_prod(3),'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
% text(model_predicted_ProLow_pu1,norm_D7_prob_ProLow,"WT PU.1 production in Pro Low",'FontSize',12,'Color','r','HorizontalAlignment','left')

% %%%%%%% together
lgd = ["","Remv. rate","","Prod. rate"];
% lgd = [lgd,"Remv. rate in pro high","",""];
legend(lgd,'location','southeast','FontSize',lg_fs)

xlabel('\it{[pu]}')
ylabel('Normalized Rate')
set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)

if contains(plot_option,'log')
    plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
    plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
else
    plot_ub_x = [0,pu1_totest_range(end)];
    plot_ub_y = [0,ylim_prod_rate];
end
xlim(plot_ub_x)
ylim(plot_ub_y)

% title(strrep( '_', ' '))
% title(final_id)
% saveas(f,fullfile('outputs_figs',strcat('_ProHigh','_',plot_option,'.tif')))
ax = gca;
ax.Box = 'on';
exportgraphics(f,fullfile(subdir_fig,strcat('_ProHigh','_',plot_option,'.tif')),'Resolution',300);

%% mac
clear lgd
% f = figure('Position',[100 100 800 600]);
% axis tight
subplot(2,2,4)
plot_option = 'linear';

% %%%%%%% generate growth curve
% tau_ProHigh = 11.6;
% tau_ProHigh_Q1 = 9.47;
% tau_ProHigh_Q3 = 15.8;

removal_range_mac = 1*pu1_totest_range;

plot(pu1_totest_range,removal_range_mac,'LineWidth',lw,'Color',removalcolor{4},'LineStyle',removalls)
hold on

%%%%%%% plot produciont rate for HIGH
cebp_cell_type_specific_conc_mac = [1];
mac_all_possibility = cebp_cell_type_specific_conc_mac;
mac_all_possibility_legend = cell(length(mac_all_possibility(:,1)),1);
mac_all_possibility_linestyle = cell(length(mac_all_possibility(:,1)),1);

for i = 1:length(mac_all_possibility(:,1))
    mac_all_possibility_legend{i} = sprintf('pro high, \\it[ce]\\rm = %.3f',mac_all_possibility(i,1));
    mac_all_possibility_linestyle{i} = linestyle_name{i};
end

TR_range_mac = zeros(size(pu1_totest_range));
for m = 1:length(mac_all_possibility(:,1))
    cebpa = mac_all_possibility(m,1);
    runx1 = 1;
    for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa,runx1,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
        
        D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        
        D7_curr_prod_rate = D7_prob;
        TR_range_mac(i) = D7_curr_prod_rate/D7_ss_prod_rate;

    end
    plot(pu1_totest_range,TR_range_mac,prohigh_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
    hold on

    cell_name = prohigh_all_possibility_legend{m};
    % lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
end




% %%%%%% find actual steady state
% [~,ind] = min(abs(TR_range_HIGH-removal_range_ProHigh));
% concise_ind = extract_intersection_index(ind,40);
% true_prohigh_pu1_ss = pu1_totest_range(concise_ind);
% true_prohigh_pu1_prod = TR_range_HIGH(concise_ind);

% %%%%%% mark itnersection
% plot(pu1_ss_conc(3),pu1_ss_prod(3),'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
% text(model_predicted_ProLow_pu1,norm_D7_prob_ProLow,"WT PU.1 production in Pro Low",'FontSize',12,'Color','r','HorizontalAlignment','left')

% %%%%%%% together
lgd = ["Remv. rate","Prod. rate"];
legend(lgd,'location','southeast','FontSize',lg_fs)

xlabel('\it{[pu]}')
ylabel('Norm. Rate')
set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)

if contains(plot_option,'log')
    plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
    plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
else
    plot_ub_x = [0,pu1_totest_range(end)];
    plot_ub_y = [0,ylim_prod_rate];
end
xlim(plot_ub_x)
ylim(plot_ub_y)

% title(strrep( '_', ' '))
% title(final_id)
% saveas(f,fullfile('outputs_figs',strcat('_ProHigh','_',plot_option,'.tif')))
ax = gca;
ax.Box = 'on';
exportgraphics(f,fullfile(subdir_fig,strcat('_mac','_',plot_option,'.tif')),'Resolution',300);

%% readable results
show_results_real_vs_pred_v3(pu1_ss_conc,pu1_ss_prod,'cebpaOnly_noPCfeedback','log',subdir_fig)