function ODE_predict_cerun_multiplePars(objFuncHandle,selec,final_fn,nuc_destab_location,CEPROD,err_threshold)
% this file is to do prediction, but across multiple good parameter sets
%% logistics
load('error_prop_0216.mat');

% subdir_fig = strcat('ODEoutputs_figs_',erase(final_fn,'.txt')); % too
% long!!!
subdir_fig = strcat('ODEoutputs_figs_multiPars','_',CEPROD); 
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
lg_fs = 16;
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

% err_threshold = 0.32;
good_mask = A(:, end) < err_threshold;
good_pars_mat = A(good_mask, 1:end-1);

%% parameter ensemble
if size(good_pars_mat, 1) > 1
    n_pars = size(good_pars_mat, 2);
    par_labels = { ...
        'G_{NUC}', 'G_{PU1 14}', 'G_{PU1 12}', 'G_{PU1 10}', 'G_{CEBP}', 'G_{RUNX}', 'G_{RNAP}', ...
        'G_{NUC PU1}', 'G_{NUC CEBP}', 'G_{NUC RUNX1}', 'G_{PU1 RNAP}', 'G_{CEBP RNAP}', 'G_{RUNX1 RNAP}', ...
        'G_{RUNX PU1}', 'G_{RUNX CEBP}', 'G_{PU1 CEBP}', 'promoter rescale', 'runx1 bcell', 'runx1 Pro', ...
        'Hill(\it[pu]\rm)_{Amin}', 'Hill(\it[pu]\rm)_{Amax}', 'Hill(\it[pu]\rm)_{N}', 'Hill(\it[pu]\rm)_{K}', ...
    };

    f=figure('Name', sprintf('Parameter Distributions (error < %.2f)',err_threshold), 'Position', [159, 283, 1677, 736]);
    hold on;
    for k = 1:n_pars
        vals = good_pars_mat(:, k);
        if numel(unique(vals)) == 1
            plot(k, vals(1), 'ro', 'MarkerSize', 6, 'LineWidth', 1, 'MarkerFaceColor', 'r');
            % text(k + 0.2, vals(1), sprintf('%.0g', vals(1)), 'FontSize', 9, 'Color', 'r', 'VerticalAlignment', 'middle');
        else
            [f_dens, xi] = ksdensity(vals);
            f_norm = f_dens / max(f_dens) * 0.4;
            fill([k - f_norm, fliplr(k + f_norm)], [xi, fliplr(xi)], [0.4 0.6 0.9], ...
                'FaceAlpha', 0.5, 'EdgeColor', 'none');
            plot(k * ones(size(vals)), vals, 'k.', 'MarkerSize', 6);
            % plot([k-0.2, k+0.2], [median(vals), median(vals)], 'r-', 'LineWidth', 2);
        end
    end
    h_constrained = plot(nan, nan, 'ro', 'MarkerSize', 6, 'LineWidth', 1, 'MarkerFaceColor', 'r');
    hold off;
    set(gca, 'XTick', 1:n_pars, 'XTickLabel', repmat({''}, 1, n_pars), 'FontSize', 16);
    tick_labels = repmat({''}, 1, n_pars);

    for k = 1:n_pars
        vals = good_pars_mat(:, k);
        if numel(unique(vals)) == 1
            tick_labels{k} = par_labels{k};
        else
            ypos = median(vals);
            if k == 19
                ytext = 5.5;
                va = 'bottom';
            elseif k == 20 || k == 23
                ytext = -0.1;
                va = 'top';
            elseif ypos < 0
                ytext = min(vals) - 0.05 * diff(ylim);
                va = 'top';
            else
                ytext = max(vals) + 0.05 * diff(ylim);
                va = 'bottom';
            end
            text(k, ytext, par_labels{k}, 'HorizontalAlignment', 'center', 'VerticalAlignment', va, ...
                'Rotation', 0, 'FontSize', 11, 'Interpreter', 'tex');
        end
    end
    legend(h_constrained, 'Fixed parameters', 'Location', 'Northeast', 'FontSize', 12);
    set(gca, 'XTickLabel', tick_labels, 'XTickLabelRotation', 45);
    ylabel('Parameter value');
    title(sprintf('Parameter distributions across %d sets with error < %.2f', size(good_pars_mat,1), err_threshold));
    box on;
end

disp('Adjust text positions in the figure, then press Enter to export.')
input('');
% saveas(f,fullfile(subdir_fig,'parametersetViolin.tif'))
exportgraphics(f,fullfile(subdir_fig,'parametersetViolin.tif'),'Resolution',600)

%% predict N times
N = size(good_pars_mat, 1);
all_pu1_ss_conc   = zeros(N, 4);
all_cebpa_ss_conc = zeros(N, 4);
all_runx1_ss_conc = zeros(N, 4);
all_pu1_ss_prod   = zeros(N, 3);
all_p_cebp_hill   = zeros(N, 4);
% all_p_auto_hill   = zeros(N, 4);
all_sim_d_final   = zeros(N, 12);
all_energy_table  = zeros(N, 16);

for le = 1:N
    [~, pu1_ss_conc, cebpa_ss_conc, runx1_ss_conc, pu1_ss_prod, p_cebp_hill,...
     real_d_final, sim_d_final, ~, energy_table] = ...
        objFuncHandle(selec, good_pars_mat(le,:), nuc_destab_location,CEPROD);

    all_pu1_ss_conc(le,:)   = pu1_ss_conc;
    all_cebpa_ss_conc(le,:) = cebpa_ss_conc;
    all_runx1_ss_conc(le,:) = runx1_ss_conc;
    all_pu1_ss_prod(le,:)   = pu1_ss_prod;
    all_p_cebp_hill(le,:)   = p_cebp_hill;
    % all_p_auto_hill(le,:)   = p_auto_hill;
    all_sim_d_final(le,:)   = sim_d_final;
    all_energy_table(le,:)  = energy_table;
end


%% predict cebpa production rate acorss multiple datasets
ceprod_good_pars = zeros(3,N); % three rows (cell types), 50+ columns (parameter set)
for le = 1:N
    for ct = 1:length(pu1_ss_conc(1:3))
        ceprod_good_pars(ct,le) = (cebpa_RegBy_pu1_Hill(all_p_cebp_hill(le,:),all_pu1_ss_conc(le,ct)))/...
            (cebpa_RegBy_pu1_Hill(all_p_cebp_hill(le,:),1));
    end
end

% plot predicted production rate for three cell types acorss parameter
% sets, versus experimental data
load('ct_TFPROD_constraint.mat','ce_Bcell_over_mac_bound','ce_pro_over_mac_bound')
f = figure('Position',[100 100 800 600]);

%%% scatter of predictions
bcell_scatter_x = zeros(1,N) + 1;
prolow_scatter_x = zeros(1,N) + 2;
prohigh_scatter_x = zeros(1,N) + 3;
scatter(bcell_scatter_x,ceprod_good_pars(1,:),'SizeData',150,'MarkerEdgeColor','b','LineWidth',lw)
hold on
scatter(prolow_scatter_x,ceprod_good_pars(2,:),'SizeData',150,'MarkerEdgeColor','b','LineWidth',lw)
hold on
scatter(prohigh_scatter_x,ceprod_good_pars(3,:),'SizeData',150,'MarkerEdgeColor','b','LineWidth',lw)

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

xlim([-0.01 3.8])
ylim([0 1])
xticks([1,2,3])
xticklabels({'Bcell','pro low','pro high'})
ylabel('C/EBPa production rate')

legend('Model predictions','','','RNA-Seq measurements','Location','northwest')
set(gca,'FontSize',fs)

ax = gca;
set(ax, 'FontSize', fs)
exportgraphics(f,fullfile(subdir_fig,'CEPRODprediction.tif'),'Resolution',600)

%% predict gene expression acorss multiple sets

load('Leddin_real_d_std.mat','real_d_std')
real_d_std = [real_d_std,5];
all_data_label = {'Pro10','Pro12','Pro14','Pro1410','Pro1412','Pro141210','PU12mut','CE14mut','PU12CE14mut','BPro14','BPro1412','RUN14mut'};

f = figure('Position',[100 100 800 600]);
h_sim = gobjects(N,1);
for le = 1:N
    h_sim(le) = plot(1:n_datap, all_sim_d_final(le,:), 'LineStyle','none','Marker','o','MarkerSize',ms,'MarkerEdgeColor',prodcolor,'LineWidth',lw);
    if le > 1; h_sim(le).HandleVisibility = 'off'; end
    hold on
end
h_real = errorbar(1:n_datap,real_d_final,real_d_std(selec)./2556.3,'LineStyle','none','Marker','o','MarkerSize',ms,'Color','k','LineWidth',lw);
set(gca, 'YScale', 'log')

xticks([1:n_datap])
xticklabels(all_data_label(selec))
ylabel('Gene expression')

ylim([0 1.3])
legend([h_sim(1), h_real],{'Model predictions','Experiment'},'Location','Northeast','FontSize',lg_fs)

set(gca,'FontSize',fs)
ax = gca;
ax.XAxis.FontSize = 16;
exportgraphics(f,fullfile(subdir_fig,'Gprediction.tif'),'Resolution',600)

% additionally plot the largeest two points
f = figure('Position',[100 100 200 400]);
plot(1:2,all_sim_d_final(:,5:6),'LineStyle','none','Marker','o','MarkerSize',ms,'MarkerEdgeColor',prodcolor,'LineWidth',lw)
hold on
errorbar(1:2,real_d_final(5:6),real_d_std(5:6)./2556.3,'LineStyle','none','Marker','o','MarkerSize',ms,'Color','k','LineWidth',lw)
set(gca, 'YScale', 'linear')

xticks([1:2])
xticklabels(all_data_label(5:6))
% ylabel('PU.1 gene expression')

ylim([0 1.2])
xlim([-0.1 2.5])
% legend('Model','Experiment','Location','Northwest','FontSize',lg_fs)

set(gca,'FontSize',fs)
exportgraphics(f,fullfile(subdir_fig,'Gprediction_141210.tif'),'Resolution',600)

%% plot multiple pu1 steady states + prod rates
show_results_real_vs_pred_v3_ensemble(all_pu1_ss_conc, all_pu1_ss_prod, [], 'log', subdir_fig)