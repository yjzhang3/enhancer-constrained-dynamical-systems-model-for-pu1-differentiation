function ODE_prediction_cebpa_runx1_PU1CEfeedback(objFuncHandle,selec,best_pars,nuc_destab_location,CEPROD,x0,subdir_fig)
%% housekeeping
load('error_prop_0216.mat');

% subdir_fig = strcat('ODEoutputs_figs_',erase(final_fn,'.txt')); % too
% long!!!
% subdir_fig = strcat('ODEoutputs_figs','_',CEPROD); 
% subdir_mat = 'ODEoutputs_mats';

% Ensure subdirectory exists (creates if not)
if ~exist(subdir_fig, 'dir')
    mkdir(subdir_fig);
end

% if ~exist(subdir_mat, 'dir')
%     mkdir(subdir_mat);
% end

RNAPFuncHandle = @RNAP_occ_prob_cebpa_runx1_general;
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
removalcolor = {'#FF5800','#BFBFBF','#808080','#2626B2'};
% removalcolor = [0.70 0.40 0.85];
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
% nvars = length(A(1,:))-1;
% sorted_arri = sort(arri);
% minerr = sorted_arri(1); % pick the smallest or second smallest?
% minind = find(arri==minerr);
% best_pars = A(minind,1:end-1);
% disp(minerr)
% previous_results = readmatrix('local_fit_with_cebpa_runx1_general_PU1CEfeedback_PrepC_constRunx1_2N_maxV1_123456789101112_010100.txt');
% [~,I] = min(previous_results(:,end));
% best_pars = previous_results(I,1:end-1);
% % best_pars = best_pars([1:17,23:26]);
% load('wRUNX1_noODEsim_ce_conc_pars_0401.mat','x0')
% x0(2,1) = 0.1734;

[final_err,pu1_ss_conc,cebpa_ss_conc,runx1_ss_conc,pu1_ss_prod,p_cebp_hill,...
    real_d_final,sim_d_final,D7_valid_sublist,energy_table] = ...
    objFuncHandle(selec,best_pars,nuc_destab_location,CEPROD,x0);
disp(final_err)
% save('Model3_ss.mat',"pu1_ss_conc","cebpa_ss_conc","runx1_ss_conc")
% save('fitting_error_by_category.mat','pu1ss_error','cebpaprod_error','pu1prod_error','gene_data_err','final_err','nvars')

%% parameter values
if ~isfile(fullfile(subdir_fig,'parametersetViolin.tif'))
    % file exists
par_labels = { ...
        'G_{NUC}', 'G_{PU1 14}', 'G_{PU1 12}', 'G_{PU1 10}', 'G_{CEBP}', 'G_{RUNX}', 'G_{RNAP}', ...
        'G_{NUC PU1}', 'G_{NUC CEBP}', 'G_{NUC RUNX1}', 'G_{PU1 RNAP}', 'G_{CEBP RNAP}', 'G_{RUNX1 RNAP}', ...
        'G_{RUNX PU1}', 'G_{RUNX CEBP}', 'G_{PU1 CEBP}', 'promoter rescale', 'runx1 bcell', 'runx1 Pro', ...
        'Hill(\it[pu]\rm)_{Amin}', 'Hill(\it[pu]\rm)_{Amax}', 'Hill(\it[pu]\rm)_{N}', 'Hill(\it[pu]\rm)_{K}', ...
    };

    f=figure('Name', 'Parameter', 'Position', [159, 283, 1677, 736]);
    hold on;
    % vals = best_pars;
    for k = 1:length(best_pars)
        % if numel(unique(vals)) == 1
            % plot(k, vals(1), 'ro', 'MarkerSize', 6, 'LineWidth', 1, 'MarkerFaceColor', 'r');
            % text(k + 0.2, vals(1), sprintf('%.0g', vals(1)), 'FontSize', 9, 'Color', 'r', 'VerticalAlignment', 'middle');
        % else
            % [f_dens, xi] = ksdensity(vals);
            % f_norm = f_dens / max(f_dens) * 0.4;
            % fill([k - f_norm, fliplr(k + f_norm)], [xi, fliplr(xi)], [0.4 0.6 0.9], ...
            %     'FaceAlpha', 0.5, 'EdgeColor', 'none');
            plot(k, best_pars(k), 'k.', 'MarkerSize', ms);
            % plot([k-0.2, k+0.2], [median(vals), median(vals)], 'r-', 'LineWidth', 2);
        % end
    end
    % h_constrained = plot(nan, nan, 'ro', 'MarkerSize', 6, 'LineWidth', 1, 'MarkerFaceColor', 'r');
    % hold off;
    set(gca, 'XTick', 1:length(best_pars), 'XTickLabel', repmat({''}, 1, length(best_pars)), 'FontSize', fs);
    tick_labels = repmat({''}, 1, length(best_pars));

    for k = 1:length(best_pars)
        % vals = good_pars_mat(:, k);
        % if numel(unique(vals)) == 1
        %     tick_labels{k} = par_labels{k};
        % else
            ypos = best_pars(k) + 0.1;
            % if k == 19
            %     ytext = 5.5;
            %     va = 'bottom';
            % elseif k == 20 || k == 23
            %     ytext = -0.1;
            %     va = 'top';
            % if ypos < 0
            %     ytext = min(best_pars) - 0.05 * diff(ylim);
            %     va = 'top';
            % else
            %     ytext = max(best_pars) + 0.05 * diff(ylim);
            %     va = 'bottom';
            % end
            text(k, ypos, par_labels{k}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'Rotation', 0, 'FontSize', lg_fs, 'Interpreter', 'tex');
        % end
    end
    % legend(h_constrained, 'Fixed parameters', 'Location', 'Northeast', 'FontSize', 12);
    set(gca, 'XTickLabel', tick_labels, 'XTickLabelRotation', 45);
    ylabel('Parameter value');
    % title(sprintf('Parameter distributions across %d sets with error < %.2f', size(good_pars_mat,1), err_threshold));
    box on;
    set(gca,'FontSize',fs,'XLim',[0 length(best_pars)+1],'YLim',[-7 7])
% end

% disp('Adjust text positions in the figure, then press Enter to export.')
% input('');
saveas(f,fullfile(subdir_fig,'parametersetViolin.tif'))
exportgraphics(f,fullfile(subdir_fig,'parametersetViolin.tif'),'Resolution',600)
end
%% plot tf Hill activation of CEBP
% f = figure('Position',[100 100 800 600]);
% 
% % compute ce production rate
% cebpa_ss_prod = zeros(1,length(pu1_ss_conc));
% % pu1_contribution = zeros(1,length(pu1_ss_conc));
% % ce_contribution = zeros(1,length(pu1_ss_conc));
% for i = 1:length(pu1_ss_conc)
%     pu1 = pu1_ss_conc(i);
%     % ce = cebpa_ss_conc(i);
%     cebpa_ss_prod(i) = cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1)/cebpa_RegBy_pu1_Hill(p_cebp_hill,1);
%     % pu1_contribution(i) = cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1)/cebpa_RegBy_pu1_Hill(p_cebp_hill,1);
%     % ce_contribution(i) = cebpa_RegBy_cebpa_Hill(p_auto_hill,ce)/cebpa_RegBy_cebpa_Hill(p_auto_hill,1);
% end
% 
% %%%%%%%%% plot ce prod rate vs. pu1 conc contribution
% 
% load('ct_TFPROD_constraint.mat','ce_Bcell_over_mac_bound','ce_pro_over_mac_bound')
% f1 = [1,2,3,4];
% 
% b1=bar([pu1_ss_conc;cebpa_ss_prod;cebpa_ss_conc;runx1_ss_conc]');
% b1(1).FaceColor = 'b';
% b1(2).FaceColor = 'm';
% b1(3).FaceColor = 'c';
% b1(4).FaceColor = 'y';
% xticks([1,2,3,4])
% xticklabels({'B cell','Pro low','Pro high','mac'})
% hold on
% pu1_ss_conc_position = b1(1).XEndPoints;
% yneg = abs(conc_ratio_all(:,1)'-conc_ratio_all(:,2)');
% ypos = abs(conc_ratio_all(:,3)'-conc_ratio_all(:,2)');
% errorbar(pu1_ss_conc_position(1:3),conc_ratio_all(:,2)',yneg,ypos,'LineStyle','none','LineWidth',lw,'Color','k')
% hold on
% 
% cebpa_ss_prod_pos = b1(2).XEndPoints;
% refpoint = 0.01;
% Bcell_yneg = abs(ce_Bcell_over_mac_bound(1)-refpoint);
% Bcell_ypos = abs(ce_Bcell_over_mac_bound(2)-refpoint);
% errorbar(cebpa_ss_prod_pos(1),refpoint,Bcell_yneg,Bcell_ypos,'LineStyle','none','LineWidth',lw,'Color','k')
% hold on
% refpoint = 0.5;
% pro_yneg = abs(ce_pro_over_mac_bound(1)-refpoint);
% pro_ypos = abs(ce_pro_over_mac_bound(2)-refpoint);
% errorbar(cebpa_ss_prod_pos(2),refpoint,pro_yneg,pro_ypos,'LineStyle','none','LineWidth',lw,'Color','k')
% hold on
% errorbar(cebpa_ss_prod_pos(3),refpoint,pro_yneg,pro_ypos,'LineStyle','none','LineWidth',lw,'Color','k')
% 
% hold on
% load('ct_TFPROD_constraint.mat','run_Bcell_over_mac_bound','run_pro_over_mac_bound')
% runx1_ss_conc_position = b1(4).XEndPoints;
% refpoint = 1;
% yneg_Bcell = abs(run_Bcell_over_mac_bound(1)-refpoint);
% ypos_Bcell = abs(run_Bcell_over_mac_bound(2)-refpoint);
% yneg_pro = abs(run_pro_over_mac_bound(1)-refpoint);
% ypos_pro = abs(run_pro_over_mac_bound(2)-refpoint);
% errorbar(runx1_ss_conc_position(1),refpoint,yneg_Bcell,ypos_Bcell,'LineStyle','none','LineWidth',lw,'Color','k')
% hold on
% errorbar(runx1_ss_conc_position(2),refpoint,yneg_pro,ypos_pro,'LineStyle','none','LineWidth',lw,'Color','k')
% hold on
% errorbar(runx1_ss_conc_position(3),refpoint,yneg_pro,ypos_pro,'LineStyle','none','LineWidth',lw,'Color','k')
% 
% legend('\it{[pu]}',strcat('Hill(\it{[pu]}) = prod(\it{[ce]})',sprintf(', N=%.3f',p_cebp_hill(3))),'\it{[ce]}','\it{[run]}','Location','northoutside')
% % ylim([0 1.2])
% ylim([10^(-2) 10^1])
% set(gca,'FontSize',fs,'YScale','log')
% saveas(f,fullfile(subdir_fig,'CEPRODprediction.tif'))

%% plot tf concentration

% f = figure('Position',[100 100 800 600]);
% 
% % conc constraint for runx1
% load('ct_TFPROD_constraint.mat','run_Bcell_over_mac_bound','run_pro_over_mac_bound')
% 
% % % range for runx1 conc in B cell
% % yline(run_Bcell_over_mac_bound(1),'LineWidth',lw,'Color',bcellcolor,'LineStyle','--')
% % hold on
% % yline(run_Bcell_over_mac_bound(2),'LineWidth',lw,'Color',bcellcolor,'LineStyle','--')
% % 
% % % range for runx1 conc in pro
% % yline(run_pro_over_mac_bound(1),'LineWidth',lw,'Color',prolowcolor,'LineStyle',':')
% % hold on
% % yline(run_pro_over_mac_bound(2),'LineWidth',lw,'Color',prolowcolor,'LineStyle',':')
% % hold on
% 
% % report conc values
% b=bar([0,1,2],[pu1_ss_conc;cebpa_ss_conc;runx1_ss_conc]);
% % colororder("reef")
% xlim([-1 3])
% 
% % Assign colors to each bar position (four bars in each group)
% b(1).FaceColor = bcellcolor;   % red
% b(2).FaceColor = prolowcolor ;% blue
% b(3).FaceColor = prohighcolor ;  % green
% b(4).FaceColor = maccolor;  % black
% b(1).FaceAlpha = fa;b(2).FaceAlpha = fa;b(3).FaceAlpha = fa;b(4).FaceAlpha = fa;
% 
% hold on
% % Loop over each bar series
% for k = 1:length(b) 
% 
%     % Endpoints of bars
%     x = b(k).XEndPoints; % x coord of the kth bar in each group of four bars
%     y = b(k).YEndPoints;
% 
%     % Corresponding values
%     labels = b(k).YData;
% 
%     % Add text labels
%     for j = 1:length(x) % jth TF
% 
%         text(x(j) - 0.05 , ...
%              y(j) + 0.05, ...
%              sprintf('%.4f', labels(j)), ...
%              'VerticalAlignment', 'middle','FontSize',17);
%         hold on
%         if (isequal(j,3) && (isequal(k,1)))
%             errorbar(x(j),runx1_ss_conc(k),abs(run_Bcell_over_mac_bound(1)-runx1_ss_conc(k)),abs(run_Bcell_over_mac_bound(2)-runx1_ss_conc(k)),'LineWidth',lw,'Color','k')
%         elseif (isequal(j,3) && (isequal(k,2))) || (isequal(j,3) && (isequal(k,3)))
%             errorbar(x(j),runx1_ss_conc(k),abs(run_pro_over_mac_bound(1)-runx1_ss_conc(k)),abs(run_pro_over_mac_bound(2)-runx1_ss_conc(k)),'LineWidth',lw,'Color','k')
%         end
%         hold on
%     end
% end
% 
% % ylabel('\it{[tf]}')
% xticks([0 1 2])
% xticklabels({'\it{[pu]}','\it{[ce]}','\it{[run]}'})
% legend('Bcell','pro low','pro high','mac','exper. [run] B cell','exper. [run] pro','FontSize',lg_fs,'Location', ...
%     'northwest')
% set(gca,'FontSize',fs)
% ax = gca;
% ax.YScale = 'log';
% 
% saveas(f,fullfile(subdir_fig,'Concprediction.tif'))

%% predict cebpa production rate acorss multiple datasets
ceprod_good_pars = zeros(length(pu1_ss_conc(1:3)),1); % three rows (cell types), 50+ columns (parameter set)
% for le = 1:N
for ct = 1:length(pu1_ss_conc(1:3))
    ceprod_good_pars(ct,1) = cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_ss_conc(ct))/(cebpa_RegBy_pu1_Hill(p_cebp_hill,1));
end
% end

% plot predicted production rate for three cell types acorss parameter
% sets, versus experimental data
load('ct_TFPROD_constraint.mat','ce_Bcell_over_mac_bound','ce_pro_over_mac_bound')
f = figure('Position',[100 100 800 600]);

%%% scatter of predictions
N = 1;
bcell_scatter_x = zeros(1,N) + 1;
prolow_scatter_x = zeros(1,N) + 2;
prohigh_scatter_x = zeros(1,N) + 3;
plot(bcell_scatter_x,ceprod_good_pars(1,:),'LineStyle','none','Marker','o','MarkerEdgeColor',prodcolor,'MarkerFaceColor',prodcolor,'LineWidth',lw,'MarkerSize',ms);
hold on
plot(prolow_scatter_x,ceprod_good_pars(2,:),'LineStyle','none','Marker','o','MarkerEdgeColor',prodcolor,'MarkerFaceColor',prodcolor,'LineWidth',lw,'MarkerSize',ms);
hold on
plot(prohigh_scatter_x,ceprod_good_pars(3,:),'LineStyle','none','Marker','o','MarkerEdgeColor',prodcolor,'MarkerFaceColor',prodcolor,'LineWidth',lw,'MarkerSize',ms);

%%% errorbar of data
hold on
refpoint = 0.01;
Bcell_yneg = abs(ce_Bcell_over_mac_bound(1)-refpoint);
Bcell_ypos = abs(ce_Bcell_over_mac_bound(2)-refpoint);
errorbar(bcell_scatter_x,refpoint,Bcell_yneg,Bcell_ypos,'LineStyle','none','LineWidth',lw,'Color','k')
hold on
refpoint = 0.5;
pro_yneg = abs(ce_pro_over_mac_bound(1)-refpoint);
pro_ypos = abs(ce_pro_over_mac_bound(2)-refpoint);
errorbar(prolow_scatter_x,refpoint,pro_yneg,pro_ypos,'LineStyle','none','LineWidth',lw,'Color','k')
hold on
errorbar(prohigh_scatter_x,refpoint,pro_yneg,pro_ypos,'LineStyle','none','LineWidth',lw,'Color','k')

% average of prediction
plot(prolow_scatter_x+0.5,(ceprod_good_pars(2,:)+ceprod_good_pars(3,:))/2,'LineWidth',lw,'LineStyle','none','Marker','pentagram','MarkerEdgeColor',prodcolor,'MarkerFaceColor',prodcolor,'MarkerSize',ms)

xlim([-0.01 3.8])
ylim([0 1])
xticks([1,2,3])
xticklabels({'Bcell','pro low','pro high'})
ylabel('C/EBPa synthesis rate')

legend('Model','','','Experiment','Location','northwest','FontSize',lg_fs)
set(gca,'FontSize',fs)

ax = gca;
% set(ax, 'FontSize', fs)
ax.Box = 'on';
exportgraphics(f,fullfile(subdir_fig,'CEPRODprediction.tif'),'Resolution',600)


%% plot prediction vs. experimental gene expression
n_datap = length(selec);
load('Leddin_real_d_std.mat','real_d_std')
real_d_std = [real_d_std,5];
all_data_label = {'Prom10','Prom12','Prom14','Prom1410','Prom1412','Prom141210','PU12mut','CE14mut','PU12CE14mut','BProm14','BProm1412','RUN14mut'};

f = figure('Position',[100 100 800 600]);
plot(1:n_datap,sim_d_final,'LineStyle','none','Marker','o','MarkerSize',ms,'MarkerFaceColor',prodcolor,'MarkerEdgeColor',prodcolor)
hold on
errorbar(1:n_datap,real_d_final,real_d_std(selec)./2556.3,'LineStyle','none','Marker','o','MarkerSize',ms,'Color','k','LineWidth',lw)
set(gca, 'YScale', 'log')

xticks([1:n_datap])
xticklabels(all_data_label(selec))
ylabel('PU.1 expression')

set(gca,'FontSize',fs)
% destab_id = num2str(nuc_destab_location);
% final_id = strrep(destab_id," ","");

% legend(strcat('Model','(',final_id,')'),'Experiment','Location','Northwest')
ylim([0 1.3])
legend('Model','Experiment','Location','Northwest','FontSize',lg_fs)

% title(strrep(orig_name, '_', ' '))
% title(final_id)
% notxt = erase(final_fn,".txt"); 
saveas(f,fullfile(subdir_fig,'Gprediction.tif'))
% close(f)
%% plot dynamics (pu vs. time) of cell-type specific model
% f = figure('Position',[100 100 800 600]);
% axis tight
% plot(tB,Bcell_over_time(:,1),'LineWidth',lw,'Color',bcellcolor)
% hold on
% plot(tB,Bcell_over_time(:,2),'LineWidth',lw,'Color',bcellcolor,'LineStyle','--')
% hold on
% plot(tPL,ProLow_over_time(:,1),'LineWidth',lw,'Color',prolowcolor)
% hold on
% plot(tPL,ProLow_over_time(:,2),'LineWidth',lw,'Color',prolowcolor,'LineStyle','--')
% hold on
% plot(tPH,ProHigh_over_time(:,1),'LineWidth',lw,'Color',prohighcolor)
% hold on
% plot(tPH,ProHigh_over_time(:,2),'LineWidth',lw,'Color',prohighcolor,'LineStyle','--')
% 
% xlabel('t')
% ylabel('\it{[pu]}')
% ylim([0 1])
% legend('B cell (pu)','B cell (ce)','pro low (pu)','pro low (ce)','pro high (pu)','pro high (ce)','FontSize',lg_fs)
% set(gca,'FontSize',fs)
% exportgraphics(f,fullfile(subdir_fig,'TimeDynamics.tif'),'Resolution',300);

%% check stability
% ODE_func_Bcell = @(t,y) ODE_func_cellType(t,y,RNAPFuncHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(1,2),remRate_constant_ratio_all_CEBPA(1,2),runx1_ss_conc(1));
% ODE_func_ProLow = @(t,y) ODE_func_cellType(t,y,RNAPFuncHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(2,2),remRate_constant_ratio_all_CEBPA(2,2),runx1_ss_conc(2));
% ODE_func_ProHigh = @(t,y) ODE_func_cellType(t,y,RNAPFuncHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,remRate_constant_ratio_all(3,1),remRate_constant_ratio_all_CEBPA(3,1),runx1_ss_conc(3));
% 
% stability_score_bc = compute_Jacobian_eig_numerical(ODE_func_Bcell,[pu1_ss_conc(1),cebpa_ss_conc(1)]);
% % disp('Bcell')
% % disp(stability_score)
% % this is to test if my function to directly compute df1dt and df2dt is
% % corect. it is!
% % fun_dF1dt(remRate_constant_ratio_all(1,2),best_pars,cebpa_ss_conc(1),pu1_ss_conc(1),runx1_ss_conc(1)) 
% % fun_dF2dt(remRate_constant_ratio_all_CEBPA(1,2),best_pars,pu1_ss_conc(1),cebpa_ss_conc(1))
% 
% stability_score_pl = compute_Jacobian_eig_numerical(ODE_func_ProLow,[pu1_ss_conc(2),cebpa_ss_conc(2)]);
% % disp('pro low')
% % disp(stability_score)
% 
% stability_score_ph = compute_Jacobian_eig_numerical(ODE_func_ProHigh,[pu1_ss_conc(3),cebpa_ss_conc(3)]);
% % disp('pro high')
% % disp(stability_score)
% 
% f = figure('Position',[100 100 800 600]);
% y = [max(stability_score_bc);max(stability_score_pl);max(stability_score_ph)];
% b=bar(y);
% % One color per group (4 rows = 4 groups)
% groupColors = [hex2rgb(bcellcolor);   
%                hex2rgb(prolowcolor);  
%                [0 0 0]]; 
% 
% nSeries = size(y, 1);
% nBarsPerSeries = size(y, 2);
% 
% for k = 1:nBarsPerSeries
%         b(k).FaceColor = 'flat';
%         for i = 1:nSeries
%             b(k).CData(i, :) = groupColors(i, :);
%         end
% end
% 
% % b(2).LineWidth = lw;       % optional: thicker edge so it's visible
% % b(2).FaceColor = 'none';
% % b(2).EdgeColor = 'flat';
% xticks([1 2 3])
% xticklabels({'B cell','pro low','pro high'})
% ylabel('Max(Re(Eigenvalues))')
% % legend('\\it[pu]\\rm','\\it[ce]\\rm','FontSize',lg_fs)
% set(gca,'FontSize',fs)
% 
% exportgraphics(f,fullfile(subdir_fig,'eigV.tif'),'Resolution',300);

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

%% all possible pu1 level and set up figure
pu1_totest_range = [0:0.0001:xlim_conc];
D7_ss_prod_rate = sim_d_final(6);
f = figure('Position',[81,64,1286,996]);
axis tight

%% bcell
% f = figure('Position',[100 100 800 600]);
% axis tight
subplot(2,2,1)
plot_option = 'log';

% %%%%%%% plot overaall color patch
patch('Faces',f1,'Vertices',color_patch_pu1_bcell_overall,'FaceColor',bcellcolor,'EdgeColor','none','FaceAlpha',fa*0.75)
hold on
lgd = [];
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
    runx1 = runx1_ss_conc(1);
    for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa,runx1,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
        
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
lgd = ["","Remv. rate","",sprintf("Prod. rate ([ce] = %.3f)",cebpa_ss_conc(1))];
legend(lgd,'location','northwest','FontSize',lg_fs)
xlabel('\it{[pu]}')
ylabel('Rate')
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
% f = figure('Position',[100 100 800 600]);
% axis tight
subplot(2,2,2)
plot_option = 'log';

% %%%%%%% plot overall patch
patch('Faces',f1,'Vertices',color_patch_pu1_ProLow_overall,'FaceColor',prolowcolor,'EdgeColor',prolowcolor,'FaceAlpha',fa)
hold on
% lgd = "Data-constrained pro low";
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
    runx1 = runx1_ss_conc(2);
    for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa,runx1,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
        
        D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        
        D7_curr_prod_rate = D7_prob;
        TR_range_LOW(i) = D7_curr_prod_rate/D7_ss_prod_rate;

    end
    plot(pu1_totest_range,TR_range_LOW,prolow_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
    hold on

    cell_name = prolow_all_possibility_legend{m};
    % lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
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
lgd = ["","Remv. rate","",sprintf("Prod. rate ([ce] = %.3f)",cebpa_ss_conc(2))];
legend(lgd,'location','northwest','FontSize',lg_fs)

xlabel('\it{[pu]}')
ylabel('Rate')
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
% f = figure('Position',[100 100 800 600]);
% axis tight
subplot(2,2,3)
plot_option = 'linear';

% %%%%%%% plot overall patch
patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',prohighcolor,'EdgeColor',prohighcolor,'FaceAlpha',fa)
hold on
% lgd = ["Data-constrained pro high"];
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
    runx1 = runx1_ss_conc(3);
    for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa,runx1,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
        
        D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        
        D7_curr_prod_rate = D7_prob;
        TR_range_HIGH(i) = D7_curr_prod_rate/D7_ss_prod_rate;

    end
    plot(pu1_totest_range,TR_range_HIGH,prohigh_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
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
lgd = ["","Remv. rate","",sprintf("Prod. rate ([ce] = %.3f)",cebpa_ss_conc(3))];
legend(lgd,'location','southeast','FontSize',lg_fs)

xlabel('\it{[pu]}')
ylabel('Rate')
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
ylabel('Rate')
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
show_results_real_vs_pred_v3(pu1_ss_conc,pu1_ss_prod,'Model3','log',subdir_fig)