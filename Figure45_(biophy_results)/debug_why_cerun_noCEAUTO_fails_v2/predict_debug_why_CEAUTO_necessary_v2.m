clear
clc
final_fn = "localfixedG_1Hill_rescaleLS1.txt";
num_Hill = 1;
selec = 1:12;
nuc_destab_location = [0,1,0,1,0,0];
subdir_fig = 'Model3_rescalLS1';

% Ensure subdirectory exists (creates if not)
if ~exist(subdir_fig, 'dir')
    mkdir(subdir_fig);
end

% if num_Hill == 1
%     fn = 'fixedG_oneHill.txt';
    objFuncHandle=@ODE_objective_cerun_0618_fixedG_find1HillPars;
        par_labels = { 'Hill(\it[pu]\rm))_{Amin}', 'Hill(\it[pu]\rm))_{Amax}', 'Hill(\it[pu]\rm))_{N}','Hill(\it[pu]\rm))_{K}'};

% elseif num_Hill == 2
%     fn = 'fixedG_twoHill.txt';
%     objFuncHandle=@ODE_objective_0618_withBestG_solvefor2HillPars;
%     par_labels = { 'Hill(\it[pu]\rm))_{Amin}', 'Hill(\it[pu]\rm))_{Amax}', 'Hill(\it[pu]\rm))_{N}','Hill(\it[pu]\rm))_{K}',...
%         'Hill(\it[ce]\rm))_{Amin}', 'Hill(\it[ce]\rm))_{Amax}', 'Hill(\it[ce]\rm))_{N}','Hill(\it[ce]\rm))_{K}'};
% end

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

[minerr,minind] = min(A(:,end));
best_pars = A(minind,1:end-1);
disp(minerr)


 [final_err,pu1_ss_conc,cebpa_ss_conc,runx1_ss_conc,pu1_ss_prod,real_d_final,sim_d_final,...
     tB,Bcell_conc_over_time,tPL,ProLow_conc_over_time,tPH,ProHigh_conc_over_time] =...
    objFuncHandle(selec,best_pars,nuc_destab_location);
disp(final_err)

%%
f=figure('position',[148,278,1088,420]);
subplot(1,2,1)
plot(tB,Bcell_conc_over_time(:,2),'LineWidth',3,'Color',bcellcolor)
hold on
plot(tPL,ProLow_conc_over_time(:,2),'LineWidth',3,'Color',prolowcolor)
hold on
plot(tPH,ProHigh_conc_over_time(:,2),'LineWidth',3,'Color',prohighcolor)
legend('Bcell','ProLow','ProHigh')
xlabel('time')
ylabel('[ce]')
set(gca,'FontSize',20)

subplot(1,2,2)
plot(tB,Bcell_conc_over_time(:,1),'LineWidth',3,'Color',bcellcolor)
hold on
plot(tPL,ProLow_conc_over_time(:,1),'LineWidth',3,'Color',prolowcolor)
hold on
plot(tPH,ProHigh_conc_over_time(:,1),'LineWidth',3,'Color',prohighcolor)
legend('Bcell','ProLow','ProHigh')
xlabel('time')
ylabel('[pu]')
set(gca,'FontSize',20)
exportgraphics(f,fullfile(subdir_fig,'fixedG_oneHill_dynamics.tif'),'Resolution',300)

%% compare G pars here to model 2
err_threshold = 10;
good_mask = A(:, end) < err_threshold;
good_pars_mat = A(good_mask, 1:end-1);

% parameter ensemble
if size(good_pars_mat, 1) > 1
    n_pars = size(good_pars_mat, 2);

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
            % if k == 19
            %     ytext = 5.5;
            %     va = 'bottom';
            % elseif k == 20 || k == 23
            %     ytext = -0.1;
            %     va = 'top';
            if ypos < 0
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

% disp('Adjust text positions in the figure, then press Enter to export.')
% input('');
% saveas(f,fullfile(subdir_fig,'parametersetViolin.tif'))
% exportgraphics(f,'parametersetViolin.tif','Resolution',600)

%% plot prediction vs. experimental gene expression
ms = 15;
load('Leddin_real_d_std.mat','real_d_std')
real_d_std = [real_d_std,10];
all_data_label = {'Prom10','Prom12','Prom14','Prom1410','Prom1412','Prom141210','PU12mut','CE14mut','PU12CE14mut','BProm14','BProm1412','RUN14mut'};

f = figure('Position',[100 100 800 600]);
plot(1:n_datap,sim_d_final,'LineStyle','none','Marker','o','MarkerSize',ms,'MarkerFaceColor',prodcolor,'MarkerEdgeColor',prodcolor)
hold on
errorbar(1:n_datap,real_d_final,real_d_std(selec)./2556.3,'LineStyle','none','Marker','o','MarkerSize',ms,'Color','k','LineWidth',lw)
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
% saveas(f,'Gprediction.tif')
exportgraphics(f,fullfile(subdir_fig,'fixedG_oneHill_Gprediction.tif'),'Resolution',300)

%% readable results
show_results_real_vs_pred_v3(pu1_ss_conc,pu1_ss_prod,'fixedG_oneHill_Gprediction','log',subdir_fig)