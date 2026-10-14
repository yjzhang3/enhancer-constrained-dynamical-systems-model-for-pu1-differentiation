clear
clc
final_fn = "noPCfeedback.txt";
err_threshold = 0.32;
selec = 1:11;
nuc_destab_location = [0,1,0,1];
subdir_name = 'noPCfeedback_results';

% Ensure subdirectory exists (creates if not)
if ~exist(subdir_name, 'dir')
    mkdir(subdir_name);
end

%% plotting housekeeping
linestyle_name = {'-','--',':','-o','->','-v','-square'};
lw = 3; % linewidth
ms = 15; % marker size
fs = 25; % fontsize
lg_fs = 18;
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

% load('wRUNX1_noODEsim_ce_conc_pars_0401.mat','x0')
% x0(2,1) = 0.1734;

%% load file%% load best pars
if ~isfile(final_fn)
    disp('file does not exist')
    good_pars = 0;
    D7_prob_bcell = 0;
    minerr = NaN;
    return
else
    n_datap = length(selec); % number of data points to fit
    A = readmatrix(final_fn);
    
end

arri = A(:,end);
sorted_arri = sort(arri);
minerr = sorted_arri(1); % pick the smallest or second smallest?
minind = find(arri==minerr);
best_pars = A(minind,1:end-1);
disp(minerr)

objFuncHandle=@ODE_objective_0618_noPCfeedback;

[final_err,noPCfeedback_pu1_ss_conc,noPCfeedback_cebpa_ss_conc,pu1_ss_prod,energy_table,...
    real_d_final,sim_d_final,D7_valid_sublist]=...
    objFuncHandle(selec,best_pars,nuc_destab_location);
disp(final_err)
% save('nofeedback_best_results.mat','best_pars','pu1_ss_conc','cebpa_ss_conc')

% %% compare G pars here to model 2
% % err_threshold = 1;
% good_mask = A(:, end) < err_threshold;
% good_pars_mat = A(good_mask, 1:end-1);
% 
% % parameter ensemble
% if size(good_pars_mat, 1) > 0
%     n_pars = size(good_pars_mat, 2);
%     par_labels = { ...
%         'G_{NUC}', 'G_{PU1 14}', 'G_{PU1 12}', 'G_{PU1 10}', 'G_{CEBP}', 'G_{RNAP}', ...
%         'G_{NUC PU1}', 'G_{NUC CEBP}', 'G_{PU1 RNAP}', 'G_{CEBP RNAP}', ...
%          'G_{PU1 CEBP}', 'promoter rescale', ...
%         '\it[ce]\rm)_{Bcell}', '\it[ce]\rm)_{Pro Low}', '\it[ce]\rm)_{Pro High}'};
% 
%     f=figure('Name', sprintf('Parameter Distributions (error < %.2f)',err_threshold), 'Position', [159, 283, 1677, 736]);
%     hold on;
%     for k = 1:n_pars
%         vals = good_pars_mat(:, k);
%         if numel(unique(vals)) == 0
%             plot(k, vals(1), 'ro', 'MarkerSize', 6, 'LineWidth', 1, 'MarkerFaceColor', 'r');
%             % text(k + 0.2, vals(1), sprintf('%.0g', vals(1)), 'FontSize', 9, 'Color', 'r', 'VerticalAlignment', 'middle');
%         else
%             [f_dens, xi] = ksdensity(vals);
%             f_norm = f_dens / max(f_dens) * 0.4;
%             fill([k - f_norm, fliplr(k + f_norm)], [xi, fliplr(xi)], [0.4 0.6 0.9], ...
%                 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%             plot(k * ones(size(vals)), vals, 'k.', 'MarkerSize', 6);
%             % plot([k-0.2, k+0.2], [median(vals), median(vals)], 'r-', 'LineWidth', 2);
%         end
%     end
%     h_constrained = plot(nan, nan, 'ro', 'MarkerSize', 6, 'LineWidth', 1, 'MarkerFaceColor', 'r');
%     hold off;
%     set(gca, 'XTick', 1:n_pars, 'XTickLabel', repmat({''}, 1, n_pars), 'FontSize', 16);
%     tick_labels = repmat({''}, 1, n_pars);
% 
%     for k = 1:n_pars
%         vals = good_pars_mat(:, k);
%         if numel(unique(vals)) == 1
%             tick_labels{k} = par_labels{k};
%         else
%             ypos = median(vals);
%             % if k == 19
%             %     ytext = 5.5;
%             %     va = 'bottom';
%             % elseif k == 20 || k == 23
%             %     ytext = -0.1;
%             %     va = 'top';
%             if ypos < 0
%                 ytext = min(vals) - 0.05 * diff(ylim);
%                 va = 'top';
%             else
%                 ytext = max(vals) + 0.05 * diff(ylim);
%                 va = 'bottom';
%             end
%             text(k, ytext, par_labels{k}, 'HorizontalAlignment', 'center', 'VerticalAlignment', va, ...
%                 'Rotation', 0, 'FontSize', 11, 'Interpreter', 'tex');
%         end
%     end
%     legend(h_constrained, 'Fixed parameters', 'Location', 'Northeast', 'FontSize', 12);
%     set(gca, 'XTickLabel', tick_labels, 'XTickLabelRotation', 45);
%     ylabel('Parameter value');
%     title(sprintf('Parameter distributions across %d sets with error < %.2f', size(good_pars_mat,1), err_threshold));
%     box on;
% end
% 
% % disp('Adjust text positions in the figure, then press Enter to export.')
% % input('');
% % saveas(f,fullfile(subdir_name,'parametersetViolin.tif'))
% exportgraphics(f,fullfile(subdir_name,'parametersetViolin.tif'),'Resolution',600)

%% plot prediction vs. experimental gene expression
ms = 15;
load('Leddin_real_d_std.mat','real_d_std')
% real_d_std = [real_d_std,10];
all_data_label = {'Prom10','Prom12','Prom14','Prom1410','Prom1412','Prom141210','PU12mut','CE14mut','PU12CE14mut','BProm14','BProm1412'};

f = figure('Position',[100 100 800 600]);
plot(1:n_datap,sim_d_final,'LineStyle','none','Marker','o','MarkerSize',ms,'MarkerFaceColor',prodcolor,'MarkerEdgeColor',prodcolor)
hold on
errorbar(1:n_datap,real_d_final,real_d_std(selec)./2556.3,'LineStyle','none','Marker','o','MarkerSize',ms,'Color','k','LineWidth',3)
set(gca, 'YScale', 'log')

xticks([1:n_datap])
xticklabels(all_data_label(selec))
ylabel('PU.1 gene expression')

set(gca,'FontSize',fs)
destab_id = num2str(nuc_destab_location);
final_id = strrep(destab_id," ","");

% legend(strcat('Model','(',final_id,')'),'Experiment','Location','Northwest')
ylim([0 1.3])
legend('Model','Experiment','Location','Northwest','FontSize',lg_fs)

% title(strrep(orig_name, '_', ' '))
% title(final_id)
notxt = erase(final_fn,".txt"); 
exportgraphics(f,fullfile(subdir_name,'noPCfeedback_Gprediction.tif'),'Resolution',600)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% plot rate vs. pu by prod and degradation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('error_prop_0216.mat')
% %% Range of [PU1] concentration and synthesis rate in which we find intersection (color patch)
% pu1_bcell_conc_range = conc_ratio_all(1,[1,3]); % skip the median of ratio because I'm plotting patch
% pu1_ProLow_conc_range = conc_ratio_all(2,[1,3]);
% pu1_ProHigh_conc_range = conc_ratio_all(3,[1,3]);
% 
% pu1_bcell_prod_rate_range = prodRate_ratio_all(1,[1,3]);
% pu1_ProLow_prod_rate_range = prodRate_ratio_all(2,[1,3]);
% pu1_ProHigh_prod_rate_range = prodRate_ratio_all(3,[1,3]);
% 
% f1 = [1 2 3 4];
% 
% % what if we combine production and [pu1] range?
% color_patch_pu1_bcell_overall = [ pu1_bcell_conc_range(1) pu1_bcell_prod_rate_range(1); pu1_bcell_conc_range(2) pu1_bcell_prod_rate_range(1); pu1_bcell_conc_range(2) pu1_bcell_prod_rate_range(2); pu1_bcell_conc_range(1) pu1_bcell_prod_rate_range(2)]; % vertices of color patch
% color_patch_pu1_ProLow_overall = [ pu1_ProLow_conc_range(1) pu1_ProLow_prod_rate_range(1); pu1_ProLow_conc_range(2) pu1_ProLow_prod_rate_range(1); pu1_ProLow_conc_range(2) pu1_ProLow_prod_rate_range(2); pu1_ProLow_conc_range(1) pu1_ProLow_prod_rate_range(2)]; % vertices of color patch
% color_patch_pu1_ProHigh_overall = [ pu1_ProHigh_conc_range(1) pu1_ProHigh_prod_rate_range(1); pu1_ProHigh_conc_range(2) pu1_ProHigh_prod_rate_range(1); pu1_ProHigh_conc_range(2) pu1_ProHigh_prod_rate_range(2); pu1_ProHigh_conc_range(1) pu1_ProHigh_prod_rate_range(2)]; % vertices of color patch
% 
% %% range of growth rate
% rem_rate_constant_Bcell = remRate_constant_ratio_all(1,2);
% rem_rate_constant_Bcell_Q1 = remRate_constant_ratio_all(1,1);
% rem_rate_constant_Bcell_Q3 = remRate_constant_ratio_all(1,3);
% 
% rem_rate_constant_ProLow = remRate_constant_ratio_all(2,2);
% rem_rate_constant_ProLow_Q1 = remRate_constant_ratio_all(2,1);
% rem_rate_constant_ProLow_Q3 = remRate_constant_ratio_all(2,3);
% 
% rem_rate_constant_ProHigh = remRate_constant_ratio_all(3,2);
% rem_rate_constant_ProHigh_Q1 = remRate_constant_ratio_all(3,1);
% rem_rate_constant_ProHigh_Q3 = remRate_constant_ratio_all(3,3);
% 
% rem_rate_constant_Mac = 1;
% 
% %% all possible pu1 level
% pu1_totest_range = [0:0.0001:xlim_conc];
% D7_ss_prod_rate = sim_d_final(6);
% f = figure('Position',[81,64,1286,996]);
% axis tight
% 
% %% bcell
% RNAPFuncHandle = @RNAP_occ_prob_cebpa_general;
% subplot(2,2,1)
% plot_option = 'log';
% 
% % %%%%%%% plot overaall color patch
% patch('Faces',f1,'Vertices',color_patch_pu1_bcell_overall,'FaceColor',bcellcolor,'EdgeColor','none','FaceAlpha',fa*0.75)
% hold on
% lgd = "Data-constrained B cell";
% hold on
% 
% % %%%%%%% generate growth curve
% % tau_Bcell = 25.5;
% % tau_Bcell_Q1 = 17.7;
% % tau_Bcell_Q3 = 41.2;
% min_removal_rate_slope = rem_rate_constant_Bcell_Q1;
% max_removal_rate_slope = rem_rate_constant_Bcell_Q3;
% 
% removal_range = rem_rate_constant_Bcell*pu1_totest_range;
% min_removal_rate = min_removal_rate_slope*pu1_totest_range;
% max_removal_rate = max_removal_rate_slope*pu1_totest_range;
% 
% plot(pu1_totest_range,removal_range,'LineWidth',lw,'Color',removalcolor{1},'LineStyle',removalls)
% hold on
% plot_range_growthrate_patch(pu1_totest_range,max_removal_rate,min_removal_rate,hex2rgb(removalcolor{1}))
% hold on
% 
% % %%%%%%% plot production rate
% cebp_cell_type_specific_conc_bcell = [cebpa_ss_conc(1)];
% bcell_all_possbility = cebp_cell_type_specific_conc_bcell;
% bcell_all_possibility_legend = cell(length(bcell_all_possbility(:,1)),1);
% bcell_all_possibility_linestyle = cell(length(bcell_all_possbility(:,1)),1);
% for i = 1:length(bcell_all_possbility(:,1))
%     bcell_all_possibility_legend{i} = sprintf('B cell, \\it[ce]\\rm = %.3f',bcell_all_possbility(i,1));
%     bcell_all_possibility_linestyle{i} = linestyle_name{i};
% end
% 
% TR_range = zeros(size(pu1_totest_range));
% 
% for m = 1:length(bcell_all_possbility(:,1))
%     cebpa = bcell_all_possbility(m,1);
%     % runx1 = bcell_all_possbility(m,2);
%     for i = 1:length(pu1_totest_range)
%         pu1 = pu1_totest_range(i);
%         totest_concentration_table = [1,pu1,cebpa,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
% 
%         D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
% 
%         D7_curr_prod_rate = D7_prob;
% 
%         TR_range(i) = D7_curr_prod_rate/D7_ss_prod_rate;
%     end
% 
%     plot(pu1_totest_range,TR_range,bcell_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
%     hold on
% 
%     cell_name = bcell_all_possibility_legend{m};
%     % lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
% end
% 
% % %%%%%% find actual steady state
% % [~,ind] = min(abs(TR_range-removal_range));
% % concise_ind = extract_intersection_index(ind,40);
% % true_Bcell_pu1_ss = pu1_totest_range(concise_ind);
% % true_Bcell_pu1_prod = TR_range(concise_ind);
% 
% % %%%%%% mark itnersection
% % xline(model_predicted_Bcell_pu1,'Color','r','LineWidth',2)
% % yline(D7_prob_bcell/D7_ss_prod_rate,'Color','r','LineWidth',2)
% % plot(pu1_ss_conc(1),pu1_ss_prod(1),'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
% % text(model_predicted_Bcell_pu1,norm_D7_prob_bcell,"WT PU.1 production in B cell",'FontSize',12,'Color','r','HorizontalAlignment','left')
% 
% % %%%%%% together
% % lgd = [lgd,"Remv. rate in B cell","",""];
% lgd = ["","Remv. rate","",sprintf("Prod. rate ([ce] = %.3f)",cebpa_ss_conc(1))];
% legend(lgd,'location','northwest','FontSize',lg_fs)
% xlabel('\it{[pu]}')
% ylabel('Rate')
% set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)
% 
% if contains(plot_option,'log')
%     plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
%     plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
% else
%     plot_ub_x = [0,pu1_totest_range(end)];
%     plot_ub_y = [0,ylim_prod_rate];
% end
% xlim(plot_ub_x)
% ylim(plot_ub_y)
% 
% % title(strrep( '_', ' '))
% % title(final_id)
% % saveas(f,fullfile(subdir_name,strcat('_Bcellonly','_',plot_option,'.tif')))
% ax = gca;
% ax.Box = 'on';
% exportgraphics(f,fullfile(subdir_name,strcat('_Bcellonly','_',plot_option,'.tif')),'Resolution',300);
% 
% 
% %% plot ProLow 
% clear lgd
% subplot(2,2,2)
% % f = figure('Position',[100 100 800 600]);
% % axis tight
% plot_option = 'log';
% 
% % %%%%%%% plot overall patch
% patch('Faces',f1,'Vertices',color_patch_pu1_ProLow_overall,'FaceColor',prolowcolor,'EdgeColor',prolowcolor,'FaceAlpha',fa)
% hold on
% lgd = "Data-constrained pro low";
% hold on
% 
% % %%%%%%% generate growth curve
% % tau_ProLow = 28.8;
% % tau_ProLow_Q1 = 19.9;
% % tau_ProLow_Q3 = 47.8;
% 
% removal_range_ProLow = rem_rate_constant_ProLow*pu1_totest_range;
% min_removal_range_ProLow = rem_rate_constant_ProLow_Q1*pu1_totest_range;
% max_removal_range_ProLow = rem_rate_constant_ProLow_Q3*pu1_totest_range;
% 
% 
% plot(pu1_totest_range,removal_range_ProLow,'LineWidth',lw,'Color',removalcolor{2},'LineStyle',removalls)
% hold on
% plot_range_growthrate_patch(pu1_totest_range,max_removal_range_ProLow,min_removal_range_ProLow,hex2rgb(removalcolor{2}))
% hold on
% 
% % %%%%%%% plot produciont rate for LOW
% cebp_cell_type_specific_conc_prolow = [cebpa_ss_conc(2)];
% prolow_all_possibility = cebp_cell_type_specific_conc_prolow;
% prolow_all_possibility_legend = cell(length(prolow_all_possibility(:,1)),1);
% prolow_all_possibility_linestyle = cell(length(prolow_all_possibility(:,1)),1);
% 
% for i = 1:length(prolow_all_possibility(:,1))
%     prolow_all_possibility_legend{i} = sprintf('pro low, \\it[ce]\\rm = %.3f',prolow_all_possibility(i,1));
%     prolow_all_possibility_linestyle{i} = linestyle_name{i};
% end
% TR_range_LOW = zeros(size(pu1_totest_range));
% 
% for m = 1:length(prolow_all_possibility(:,1))
%     cebpa = prolow_all_possibility(m,1);
%     % runx1 = prolow_all_possibility(m,2);
%     for i = 1:length(pu1_totest_range)
%         pu1 = pu1_totest_range(i);
%         totest_concentration_table = [1,pu1,cebpa,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
% 
%         D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
% 
%         D7_curr_prod_rate = D7_prob;
%         TR_range_LOW(i) = D7_curr_prod_rate/D7_ss_prod_rate;
% 
%     end
%     plot(pu1_totest_range,TR_range_LOW,prolow_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
%     hold on
% 
%     cell_name = prolow_all_possibility_legend{m};
%     lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
% end
% 
% 
% % %%%%%% find actual steady state
% % [~,ind] = mink(abs(TR_range_LOW-removal_range_ProLow),1000);
% % % concise_ind = extract_intersection_index(ind,30);
% % true_prolow_pu1_ss = pu1_totest_range(ind(end));
% % true_prolow_pu1_prod = TR_range_LOW(ind(end));
% 
% % %%%%%% mark itnersection
% % plot(pu1_ss_conc(2),pu1_ss_prod(2),'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
% % text(model_predicted_ProLow_pu1,norm_D7_prob_ProLow,"WT PU.1 production in Pro Low",'FontSize',12,'Color','r','HorizontalAlignment','left')
% 
% % %%%%%%% together
% % lgd = [lgd,"Remv. rate in pro low","",""];
% lgd = ["","Remv. rate","",sprintf("Prod. rate ([ce] = %.3f)",cebpa_ss_conc(2))];
% legend(lgd,'location','northwest','FontSize',lg_fs)
% 
% xlabel('\it{[pu]}')
% ylabel('Rate')
% set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)
% 
% if contains(plot_option,'log')
%     plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
%     plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
% else
%     plot_ub_x = [0,pu1_totest_range(end)];
%     plot_ub_y = [0,ylim_prod_rate];
% end
% xlim(plot_ub_x)
% ylim(plot_ub_y)
% 
% % title(strrep( '_', ' '))
% % title(final_id)
% % saveas(f,fullfile('outputs_figs',strcat('_ProLow','_',plot_option,'.tif')))
% ax = gca;
% ax.Box = 'on';
% exportgraphics(f,fullfile(subdir_name,strcat('_ProLow','_',plot_option,'.tif')),'Resolution',300);
% 
% %% then ProHigh
% clear lgd
% subplot(2,2,3)
% % f = figure('Position',[100 100 800 600]);
% % axis tight
% plot_option = 'linear';
% 
% % %%%%%%% plot overall patch
% patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',prohighcolor,'EdgeColor',prohighcolor,'FaceAlpha',fa)
% hold on
% lgd = ["Data-constrained pro high"];
% hold on
% 
% % %%%%%%% generate growth curve
% % tau_ProHigh = 11.6;
% % tau_ProHigh_Q1 = 9.47;
% % tau_ProHigh_Q3 = 15.8;
% 
% removal_range_ProHigh = rem_rate_constant_ProHigh*pu1_totest_range;
% min_removal_range_ProHigh = rem_rate_constant_ProHigh_Q1*pu1_totest_range;
% max_removal_range_ProHigh = rem_rate_constant_ProHigh_Q3*pu1_totest_range;
% 
% plot(pu1_totest_range,min_removal_range_ProHigh,'LineWidth',lw,'Color',removalcolor{3},'LineStyle',removalls)
% hold on
% plot_range_growthrate_patch(pu1_totest_range,max_removal_range_ProHigh,min_removal_range_ProHigh,hex2rgb(removalcolor{3}))
% 
% 
% %%%%%%% plot produciont rate for HIGH
% cebp_cell_type_specific_conc_prohigh = [cebpa_ss_conc(3)];
% prohigh_all_possibility = cebp_cell_type_specific_conc_prohigh;
% prohigh_all_possibility_legend = cell(length(prohigh_all_possibility(:,1)),1);
% prohigh_all_possibility_linestyle = cell(length(prohigh_all_possibility(:,1)),1);
% 
% for i = 1:length(prohigh_all_possibility(:,1))
%     prohigh_all_possibility_legend{i} = sprintf('pro high, \\it[ce]\\rm = %.3f',prohigh_all_possibility(i,1));
%     prohigh_all_possibility_linestyle{i} = linestyle_name{i};
% end
% 
% TR_range_HIGH = zeros(size(pu1_totest_range));
% for m = 1:length(prohigh_all_possibility(:,1))
%     cebpa = prohigh_all_possibility(m,1);
%     % runx1 = prohigh_all_possibility(m,2);
%     for i = 1:length(pu1_totest_range)
%         pu1 = pu1_totest_range(i);
%         totest_concentration_table = [1,pu1,cebpa,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
% 
%         D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
% 
%         D7_curr_prod_rate = D7_prob;
%         TR_range_HIGH(i) = D7_curr_prod_rate/D7_ss_prod_rate;
% 
%     end
%     plot(pu1_totest_range,TR_range_HIGH,prohigh_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
%     hold on
% 
%     cell_name = prohigh_all_possibility_legend{m};
%     lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
% end
% 
% 
% 
% % %%%%%% find actual steady state
% % [~,ind] = min(abs(TR_range_HIGH-removal_range_ProHigh));
% % concise_ind = extract_intersection_index(ind,40);
% % true_prohigh_pu1_ss = pu1_totest_range(concise_ind);
% % true_prohigh_pu1_prod = TR_range_HIGH(concise_ind);
% 
% % %%%%%% mark itnersection
% % plot(pu1_ss_conc(3),pu1_ss_prod(3),'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
% % text(model_predicted_ProLow_pu1,norm_D7_prob_ProLow,"WT PU.1 production in Pro Low",'FontSize',12,'Color','r','HorizontalAlignment','left')
% 
% % %%%%%%% together
% lgd = ["","Remv. rate","",sprintf("Prod. rate ([ce] = %.3f)",cebpa_ss_conc(3))];
% % lgd = [lgd,"Remv. rate in pro high","",""];
% legend(lgd,'location','southeast','FontSize',lg_fs)
% 
% xlabel('\it{[pu]}')
% ylabel('Rate')
% set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)
% 
% if contains(plot_option,'log')
%     plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
%     plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
% else
%     plot_ub_x = [0,pu1_totest_range(end)];
%     plot_ub_y = [0,ylim_prod_rate];
% end
% xlim(plot_ub_x)
% ylim(plot_ub_y)
% 
% % title(strrep( '_', ' '))
% % title(final_id)
% % saveas(f,fullfile('outputs_figs',strcat('_ProHigh','_',plot_option,'.tif')))
% ax = gca;
% ax.Box = 'on';
% exportgraphics(f,fullfile(subdir_name,strcat('_ProHigh','_',plot_option,'.tif')),'Resolution',300);
% 
% %% mac
% clear lgd
% % f = figure('Position',[100 100 800 600]);
% % axis tight
% subplot(2,2,4)
% plot_option = 'linear';
% 
% % %%%%%%% generate growth curve
% % tau_ProHigh = 11.6;
% % tau_ProHigh_Q1 = 9.47;
% % tau_ProHigh_Q3 = 15.8;
% 
% removal_range_mac = 1*pu1_totest_range;
% 
% plot(pu1_totest_range,removal_range_mac,'LineWidth',lw,'Color',removalcolor{4},'LineStyle',removalls)
% hold on
% 
% %%%%%%% plot produciont rate for HIGH
% cebp_cell_type_specific_conc_mac = [1];
% mac_all_possibility = cebp_cell_type_specific_conc_mac;
% mac_all_possibility_legend = cell(length(mac_all_possibility(:,1)),1);
% mac_all_possibility_linestyle = cell(length(mac_all_possibility(:,1)),1);
% 
% for i = 1:length(mac_all_possibility(:,1))
%     mac_all_possibility_legend{i} = sprintf('pro high, \\it[ce]\\rm = %.3f',mac_all_possibility(i,1));
%     mac_all_possibility_linestyle{i} = linestyle_name{i};
% end
% 
% TR_range_mac = zeros(size(pu1_totest_range));
% for m = 1:length(mac_all_possibility(:,1))
%     cebpa = mac_all_possibility(m,1);
%     runx1 = 1;
%     for i = 1:length(pu1_totest_range)
%         pu1 = pu1_totest_range(i);
%         totest_concentration_table = [1,pu1,cebpa,runx1,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
% 
%         D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
% 
%         D7_curr_prod_rate = D7_prob;
%         TR_range_mac(i) = D7_curr_prod_rate/D7_ss_prod_rate;
% 
%     end
%     plot(pu1_totest_range,TR_range_mac,prohigh_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
%     hold on
% 
%     cell_name = prohigh_all_possibility_legend{m};
%     % lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
% end
% 
% 
% 
% 
% % %%%%%% find actual steady state
% % [~,ind] = min(abs(TR_range_HIGH-removal_range_ProHigh));
% % concise_ind = extract_intersection_index(ind,40);
% % true_prohigh_pu1_ss = pu1_totest_range(concise_ind);
% % true_prohigh_pu1_prod = TR_range_HIGH(concise_ind);
% 
% % %%%%%% mark itnersection
% % plot(pu1_ss_conc(3),pu1_ss_prod(3),'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
% % text(model_predicted_ProLow_pu1,norm_D7_prob_ProLow,"WT PU.1 production in Pro Low",'FontSize',12,'Color','r','HorizontalAlignment','left')
% 
% % %%%%%%% together
% lgd = ["Remv. rate","Prod. rate"];
% legend(lgd,'location','southeast','FontSize',lg_fs)
% 
% xlabel('\it{[pu]}')
% ylabel('Rate')
% set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)
% 
% if contains(plot_option,'log')
%     plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
%     plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
% else
%     plot_ub_x = [0,pu1_totest_range(end)];
%     plot_ub_y = [0,ylim_prod_rate];
% end
% xlim(plot_ub_x)
% ylim(plot_ub_y)
% 
% % title(strrep( '_', ' '))
% % title(final_id)
% % saveas(f,fullfile('outputs_figs',strcat('_ProHigh','_',plot_option,'.tif')))
% ax = gca;
% ax.Box = 'on';
% exportgraphics(f,fullfile(subdir_name,strcat('_mac','_',plot_option,'.tif')),'Resolution',300);


%% readable results
show_results_real_vs_pred_v3(noPCfeedback_pu1_ss_conc,pu1_ss_prod,'noPCfeedback','log',subdir_name)

%% now find hill mapping between cebpa_ss_conc and pu1_ss_conc assuming Hill(pu)
PCn = 10;

% find a set of Hill pars that satisfy best pu1_ss map to best ce_ss
objFunc = @(Hill_pars) find_best_pars(noPCfeedback_pu1_ss_conc,noPCfeedback_cebpa_ss_conc,Hill_pars);

% use the same bound for model 3
lb = [0,1,PCn,0.00001];
ub = [0.2,1,PCn,2];
fn = strcat('OneHillPars_fitted_from_noPCfeedback',num2str(PCn),'.txt');

% for i = 1:10
% [x,fval] = particleswarm(objFunc,4,lb,ub);
% create_parameter_file(fn,x,fval,'a+')
% end

% extract best pars
A = readmatrix(fn);
[minerr,I] = min(A(:,end));
best_pars_oneHill = A(I,1:end-1);
[err,pred_ce_conc] = find_best_pars(noPCfeedback_pu1_ss_conc,noPCfeedback_cebpa_ss_conc,best_pars_oneHill);

% plot
close all
f = figure('Position',[100 100 800 600]);
plot(noPCfeedback_pu1_ss_conc,noPCfeedback_cebpa_ss_conc,'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',15)
hold on
plot(noPCfeedback_pu1_ss_conc,pred_ce_conc,'LineStyle','none','Marker','o','MarkerFaceColor','none','MarkerEdgeColor','k','MarkerSize',15,'LineWidth',lw)
legend('Results without feedback',['Results mapping \it{[ce]}\rm to \it{[pu]}\rm,',' ',sprintf('N = %.0f',PCn)],'location','SouthEast','FontSize',lg_fs)
xlabel('Fitted \it{[pu]}')
ylabel('Fitted \it{[ce]}')
set(gca,'FontSize',25)
exportgraphics(f,fullfile(subdir_name,strcat('OneHillPars_fitted_from_noPCfeedback',num2str(PCn),'.tif')),'Resolution',300)

%% only if the increase in cebpa is captured well, then now plug one hill pars 
% back into the ODE to see what happens to the stability of the 
% stable ss found by nofeedback
load('error_prop_0216.mat')
close all
% color for constraint patch
bcellcolor = '#FF5800';
prolowcolor = '#BFBFBF';
prohighcolor = 'k';
maccolor = '#2626B2';
color_scheme = {bcellcolor,prolowcolor,prohighcolor};

totest_pars = [best_pars(1:12),best_pars_oneHill];
x0_pu1_ss_conc = noPCfeedback_pu1_ss_conc;
x0_cebpa_ss_conc = noPCfeedback_cebpa_ss_conc;
[Bcell_conc_over_time,ProLow_conc_over_time,ProHigh_conc_over_time,tB,tPL,tPH,...
    real_d_final,sim_d_final,plug_PCfeedback_pu1_ss_conc,plug_PCfeedback_pu1_ss_prod]= ...
    ODE_objective_plug_PCnoAUTOfeedback_0618(totest_pars,nuc_destab_location,x0_pu1_ss_conc,x0_cebpa_ss_conc);

%% predict first
show_results_real_vs_pred_v3(plug_PCfeedback_pu1_ss_conc,plug_PCfeedback_pu1_ss_prod,'plug_PCnoAUTOfeedback','log',subdir_name)

ms = 15;
load('Leddin_real_d_std.mat','real_d_std')
% real_d_std = [real_d_std,10];
all_data_label = {'Prom10','Prom12','Prom14','Prom1410','Prom1412','Prom141210','PU12mut','CE14mut','PU12CE14mut','BProm14','BProm1412'};

f = figure('Position',[100 100 800 600]);
plot(1:n_datap,sim_d_final,'LineStyle','none','Marker','o','MarkerSize',ms,'MarkerFaceColor',prodcolor,'MarkerEdgeColor',prodcolor)
hold on
errorbar(1:n_datap,real_d_final,real_d_std(selec)./2556.3,'LineStyle','none','Marker','o','MarkerSize',ms,'Color','k','LineWidth',3)
set(gca, 'YScale', 'log')

xticks([1:n_datap])
xticklabels(all_data_label(selec))
ylabel('PU.1 gene expression')

set(gca,'FontSize',fs)
destab_id = num2str(nuc_destab_location);
final_id = strrep(destab_id," ","");

% legend(strcat('Model','(',final_id,')'),'Experiment','Location','Northwest')
ylim([0 1.3])
legend('Model','Experiment','Location','Northwest','FontSize',lg_fs)

% title(strrep(orig_name, '_', ' '))
% title(final_id)
notxt = erase(final_fn,".txt"); 
exportgraphics(f,fullfile(subdir_name,'plug_PCnoAUTOfeedback_Gprediction.tif'),'Resolution',600)

%% examine bistability
close all;lst = '--';
% f=figure('position',[148,59,1327,1057]);

f = figure('Position',[100 100 800 600]);
% subplot(2,2,1)
xl = [0 200];  % get the current x-axis limits [xmin xmax] from the open figure
for ct = 1:3
y0 = conc_ratio_all(ct,1);      % baseline y-value where the patch starts (adjust to your needs)
w = conc_ratio_all(ct,3) - conc_ratio_all(ct,1);
p = patch([xl(1) xl(2) xl(2) xl(1) ],[y0 y0 y0+w y0+w], [0.8 0.8 0.8], ...
    'FaceAlpha', fa*0.8, 'FaceColor',color_scheme{ct},'EdgeColor','none');
hold on
end

plot(tB,Bcell_conc_over_time(:,1),'Color',bcellcolor,'LineWidth',3,'LineStyle',lst)
hold on
plot(tPL,ProLow_conc_over_time(:,1),'Color',prolowcolor,'LineWidth',3,'LineStyle',lst)
hold on
plot(tPH,ProHigh_conc_over_time(:,1),'Color',prohighcolor,'LineWidth',3,'LineStyle',lst)
legend('','','','B cell','Pro Low (IC1)','Pro High','FontSize',lg_fs,'Location','best')
ylim([0 1])
xlabel('Time (hr)')
ylabel('\it[pu]')
set(gca,'FontSize',fs,'YScale','linear')
box on;
exportgraphics(f,fullfile(subdir_name,strcat('IC1_PU_Plug_OneHillPars_fitted_from_noPCfeedback',num2str(PCn),'.tif')),'Resolution',300)

f = figure('Position',[100 100 800 600]);
plot(tB,Bcell_conc_over_time(:,2),'Color',bcellcolor,'LineWidth',3,'LineStyle',lst)
hold on
plot(tPL,ProLow_conc_over_time(:,2),'Color',prolowcolor,'LineWidth',3,'LineStyle',lst)
hold on
plot(tPH,ProHigh_conc_over_time(:,2),'Color',prohighcolor,'LineWidth',3,'LineStyle',lst)
% legend('B cell','Pro Low (IC1)','Pro High','FontSize',lg_fs,'Location','best')
xlabel('Time (hr)')
ylabel('\it[ce]')
ylim([0 1])
set(gca,'FontSize',fs,'YScale','linear')
box on;
exportgraphics(f,fullfile(subdir_name,strcat('IC1_CE_Plug_OneHillPars_fitted_from_noPCfeedback',num2str(PCn),'.tif')),'Resolution',300)

%%%%%%%%%%%%%%%%%%%% now use a different initial condition
f = figure('Position',[100 100 800 600]);
totest_pars = [best_pars(1:12),best_pars_oneHill];
x0_pu1_ss_conc(2) = 0.1;
x0_cebpa_ss_conc(2) = 0.1;
[Bcell_conc_over_time2,ProLow_conc_over_time2,ProHigh_conc_over_time2,tB,tPL,tPH]= ...
    ODE_objective_plug_PCnoAUTOfeedback_0618(totest_pars,nuc_destab_location,x0_pu1_ss_conc,x0_cebpa_ss_conc);
% subplot(2,2,2)
xl = [0 200];  % get the current x-axis limits [xmin xmax] from the open figure
for ct = 1:3
y0 = conc_ratio_all(ct,1);      % baseline y-value where the patch starts (adjust to your needs)
w = conc_ratio_all(ct,3) - conc_ratio_all(ct,1);
p = patch([xl(1) xl(2) xl(2) xl(1) ],[y0 y0 y0+w y0+w], [0.8 0.8 0.8], ...
    'FaceAlpha', fa*0.8, 'FaceColor',color_scheme{ct},'EdgeColor','none');
hold on
end
plot(tB,Bcell_conc_over_time2(:,1),'Color',bcellcolor,'LineWidth',3,'LineStyle',lst)
hold on
plot(tPL,ProLow_conc_over_time2(:,1),'Color',prolowcolor,'LineWidth',3,'LineStyle',lst)
hold on
plot(tPH,ProHigh_conc_over_time2(:,1),'Color',prohighcolor,'LineWidth',3,'LineStyle',lst)
legend('','','','B cell','Pro Low (IC2)','Pro High','FontSize',lg_fs,'Location','best')
xlabel('Time (hr)')
ylabel('\it[pu]')
ylim([0 1])
set(gca,'FontSize',fs,'YScale','linear')
box on;
exportgraphics(f,fullfile(subdir_name,strcat('IC2_PU_Plug_OneHillPars_fitted_from_noPCfeedback',num2str(PCn),'.tif')),'Resolution',300)

f = figure('Position',[100 100 800 600]);
% subplot(2,2,4)
plot(tB,Bcell_conc_over_time2(:,2),'Color',bcellcolor,'LineWidth',3,'LineStyle',lst)
hold on
plot(tPL,ProLow_conc_over_time2(:,2),'Color',prolowcolor,'LineWidth',3,'LineStyle',lst)
hold on
plot(tPH,ProHigh_conc_over_time2(:,2),'Color',prohighcolor,'LineWidth',3,'LineStyle',lst)
% legend('B cell','Pro Low (IC2)','Pro High','FontSize',lg_fs,'Location','best')
xlabel('Time (hr)')
ylabel('\it[ce]')
ylim([1e-2 1])
set(gca,'FontSize',fs,'YScale','linear')
box on;
exportgraphics(f,fullfile(subdir_name,strcat('IC2_CE_Plug_OneHillPars_fitted_from_noPCfeedback',num2str(PCn),'.tif')),'Resolution',300)

% %% now find hill mapping between cebpa_ss_conc and pu1_ss_conc, PC feedback and autoregulation
% PCn = 2;
% AUTOn = 4;
% % find a set of Hill pars that satisfy best pu1_ss map to best ce_ss
% objFunc = @(p) find_best_pars_twoHill(pu1_ss_conc,cebpa_ss_conc,p);
% 
% % use the same bound for model 3
% lb = [0.0001,1,PCn,0.001,0.0001,1,AUTOn,0.001];
% ub = [0.1,1,PCn,3,0.1,1,AUTOn,3];
% 
% fn = strcat('TwoHillPars_fitted_from_noPCfeedback',num2str(PCn),num2str(AUTOn),'.txt');
% 
% for i = 1:10
% [x,fval] = particleswarm(objFunc,8,lb,ub);
% create_parameter_file(fn,x,fval,'a+')
% end
% 
% % extract best pars
% A = readmatrix(fn);
% [minerr,I] = min(A(:,end));
% best_pars_twoHill = A(I,1:end-1);
% [err,pred_ce_conc] = find_best_pars_twoHill(pu1_ss_conc,cebpa_ss_conc,best_pars_twoHill);
% 
% % plot
% f=figure();
% plot(pu1_ss_conc,cebpa_ss_conc,'LineStyle','none','Marker','o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',15)
% hold on
% plot(pu1_ss_conc,pred_ce_conc,'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',15)
% legend('Best-fitted [ce]',strcat('Predicted [ce],',sprintf('PCn=%.3f, AUTOn=%.3f',PCn,AUTOn)),'location','best')
% xlabel('Best-fitted [pu]')
% ylabel('[ce]')
% set(gca,'FontSize',20)
% exportgraphics(f,strcat('TwoHillPars_fitted_from_noPCfeedback',num2str(PCn),num2str(AUTOn),'.tif'),'resolution',300)
% 
% 
% %% only if the increase in cebpa is capture well, then plug back into the ODE
% totest_pars = [best_pars(1:12),best_pars_twoHill];
% x0_pu1_ss_conc = pu1_ss_conc;
% x0_cebpa_ss_conc = cebpa_ss_conc;
% [new_pu1_ss_conc,new_cebpa_ss_conc,tB,Bcell_conc_over_time,tPL,ProLow_conc_over_time,tPH,ProHigh_conc_over_time]= ...
%     ODE_objective_2D_cebpa_PCAUTOfeedback_smallestPHdilute_0618(totest_pars,nuc_destab_location,x0_pu1_ss_conc,x0_cebpa_ss_conc);
% 
% close all
% f=figure('position',[148,278,1088,420]);
% subplot(1,2,1)
% plot(tB,Bcell_conc_over_time(:,1),'Color',bcellcolor,'LineWidth',3)
% hold on
% plot(tPL,ProLow_conc_over_time(:,1),'Color',prolowcolor,'LineWidth',3)
% hold on
% plot(tPH,ProHigh_conc_over_time(:,1),'Color',prohighcolor,'LineWidth',3)
% legend('B cell','Pro Low','Pro High')
% xlabel('time')
% ylabel('[pu]')
% set(gca,'FontSize',20)
% 
% subplot(1,2,2)
% plot(tB,Bcell_conc_over_time(:,2),'Color',bcellcolor,'LineWidth',3)
% hold on
% plot(tPL,ProLow_conc_over_time(:,2),'Color',prolowcolor,'LineWidth',3)
% hold on
% plot(tPH,ProHigh_conc_over_time(:,2),'Color',prohighcolor,'LineWidth',3)
% legend('B cell','Pro Low','Pro High')
% xlabel('time')
% ylabel('[ce]')
% set(gca,'FontSize',20)
% exportgraphics(f,strcat('Plug_TwoHillPars_fitted_from_noPCfeedback',num2str(PCn),num2str(AUTOn),'.tif'),'resolution',300)
% 
