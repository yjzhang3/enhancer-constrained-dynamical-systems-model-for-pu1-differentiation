function show_results_real_vs_pred_v3(pred_pu1_SS,pred_pu1_prod,figurename,axisScale,subdir_fig)

%% housekeeping
% pred_pu1_SS is a 1by3 array (from left to right, B cell, pro low, pro
% high)
% pred_pu1_prod is a 1by 3 array

% SS_color = '#E03B3B';
% prod_color = [0.00 0.75 0.65];
fs = 25;
lg_fs = 18;
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
for i = 1:3
    % Extract bounds
    xMin = matA(i,1);
    xMax = matA(i,3);
    yMin = matB(i,1);
    yMax = matB(i,3);

    % Define rectangle vertices (clockwise)
    xVertices = [xMin xMax xMax xMin];
    yVertices = [yMin yMin yMax yMax];

    % Plot patch
    patch('XData',xVertices, 'YData',yVertices, 'FaceColor',patch_colors{i}, ...
          'FaceAlpha', fa, ...
          'EdgeColor', 'none', ...
          'LineWidth', lw);
end
hold on

for i = 1:3
    % plot predicted data
    plot(pred_pu1_SS(i),pred_pu1_prod(i),'LineStyle','none','Marker',pred_marker_type,'MarkerFaceColor',dot_colors{i},'MarkerEdgeColor','none','MarkerSize',ms)
    hold on
end

% plot(1,1,'LineStyle','none','Marker','pentagram','MarkerFaceColor',patch_colors{4},'MarkerEdgeColor','none','MarkerSize',ms*1.5);

% legend({'Real B cell','Real pro low','Real pro high','Predicted B cell','Predicted pro low','Predicted pro high','Macrophage (ref)'}, ...
%        'Location', 'NorthWest','FontSize',lg_fs); 
legend({'','','','B cell','Pro low','Pro high'}, ...
       'Location', 'NorthWest','FontSize',lg_fs); 

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
exportgraphics(f,fullfile(subdir_fig,strcat(figurename,'_combinedResults','.tif')),'Resolution',300);

