function show_results_real_vs_pred_v3_ensemble(all_pu1_ss_conc, all_pu1_ss_prod, figurename, axisScale, subdir_fig)

%% housekeeping
% all_pu1_ss_conc is N x 4; columns 1-3 = B cell, pro low, pro high
% all_pu1_ss_prod is N x 3; columns 1-3 = B cell, pro low, pro high

fs = 25;
lg_fs = 16;
fa = 0.6;
lw = 3;
ms = 13;
pred_marker_type = 'o';
ss_lim = 1.2;
prod_lim = 2.5;
load('error_prop_0216.mat');

matA = conc_ratio_all;
matB = prodRate_ratio_all;

patch_colors = {'#FF5800','#BFBFBF','k','#2626B2'};
dot_colors = {'#FF5800','#808080','k','#2626B2'};

%% plot
f = figure('Position',[100 100 800 600]);
axis tight
h_patch = gobjects(1, 3);
for i = 1:3
    xMin = matA(i,1);
    xMax = matA(i,3);
    yMin = matB(i,1);
    yMax = matB(i,3);
    xVertices = [xMin xMax xMax xMin];
    yVertices = [yMin yMin yMax yMax];
    h_patch(i) = patch('XData',xVertices, 'YData',yVertices, 'FaceColor',patch_colors{i}, ...
          'FaceAlpha', fa, 'EdgeColor', 'none', 'LineWidth', lw);
end
hold on

N = size(all_pu1_ss_conc, 1);
h_pred = gobjects(1, 3);
for le = 1:N
    for i = 1:3
        h = plot(all_pu1_ss_conc(le,i), all_pu1_ss_prod(le,i), ...
            'LineStyle', 'none', 'Marker', pred_marker_type, ...
            'MarkerFaceColor', 'none', 'MarkerEdgeColor', dot_colors{i}, 'MarkerSize', ms, ...
            'LineWidth', lw);
        if le == 1
            h_pred(i) = h;
        else
            h.HandleVisibility = 'off';
        end
        hold on
    end
end

h_mac = plot(1, 1, 'LineStyle', 'none', 'Marker', 'pentagram', ...
    'MarkerFaceColor', patch_colors{4}, 'MarkerEdgeColor', 'none', 'MarkerSize', ms*1.5);

legend([h_patch(1), h_patch(2), h_patch(3), h_pred(1), h_pred(2), h_pred(3), h_mac], ...
       {'Real B cell', 'Real pro low', 'Real pro high', 'Predicted B cell', 'Predicted pro low', 'Predicted pro high', 'Macrophage (ref)'}, ...
       'Location', 'NorthWest', 'FontSize', lg_fs);

%% Plot formatting
xlim([0 ss_lim])
ylim([0 prod_lim])
xlabel('\it{[pu]}');
ylabel('Synthesis rate');
ax = gca;
ax.Box = 'on';
ax.FontSize = fs;
ax.XScale = axisScale;
ax.YScale = axisScale;
exportgraphics(f, fullfile(subdir_fig, strcat(figurename, '_ensembleResults', '.tif')), 'Resolution', 300);
