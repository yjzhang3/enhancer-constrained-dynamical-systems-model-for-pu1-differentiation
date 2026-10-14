function show_results_real_vs_pred(pred_pu1_SS,pred_pu1_prod,figurename)

% pred_pu1_SS is a 1by3 array (from left to right, B cell, pro low, pro
% high)
% pred_pu1_prod is a 1by 3 array

SS_color = '#E03B3B';
prod_color = [0.00 0.75 0.65];
fs = 19;
lw = 3;
ms = 10;
load(fullfile(fileparts(pwd), 'Error propagation and Normalization', 'error_prop.mat'));

%%%% categrocial x
cats = categorical({'B cell','pro low','pro high'},'Ordinal',true);
% categories(cats) % debug
desire_order = {'B cell','pro low','pro high'};
desire_cats = reordercats(cats,desire_order);
% categories(desire_cats) % debug
X = double(desire_cats);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
f = figure('Position',[100 100 800 600]);
axis tight

%%%% plot predicted pu1 ss
bar(desire_cats,pred_pu1_SS,'FaceColor','none', 'EdgeColor',SS_color,'LineWidth',lw)
hold on

%%%% overlay predicted pu1 ss
for ct = 1:3
    curr_conc_range = conc_ratio_all(ct,:);
    errorbar(X(ct),curr_conc_range(2), ...% x, y
        curr_conc_range(2)-curr_conc_range(1),curr_conc_range(3)-curr_conc_range(2), ... % yneg, ypos
        'LineWidth',lw,'Color',SS_color,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',SS_color,'MarkerSize',ms);
    hold on
end
legend('Model-predicted \it{[pu]}','Data-constrained PU.1-GFP', ...
    'Location','northwest')

%%%% save file name
set(gca,'FontSize',fs)
exportgraphics(f,strcat(figurename,'pu1_ss','.tif'),'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure('Position',[100 100 800 600]);
axis tight

%%%% plot predicted pu1 prod
bar(desire_cats,pred_pu1_prod,'FaceColor','none', 'EdgeColor',prod_color,'LineWidth',lw)
hold on

%%%% overlay predicted pu1 prod
for ct = 1:3
    curr_prod_rate_range = prodRate_ratio_all(ct,:);
    errorbar(X(ct),curr_prod_rate_range(2), ...% x, y
        curr_prod_rate_range(2)-curr_prod_rate_range(1),curr_prod_rate_range(3)-curr_prod_rate_range(2), ... % yneg, ypos
        'LineWidth',lw,'Color',prod_color,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',prod_color,'MarkerSize',ms);
    hold on
end

legend('Model-predicted \it[pu]\rm production rate','Data-constrained PU.1 synthesis rate', ...
    'Location','northwest')

%%%% save file name
set(gca,'FontSize',fs)
exportgraphics(f,strcat(figurename,'pu1_prod','.tif'),'Resolution',300);
