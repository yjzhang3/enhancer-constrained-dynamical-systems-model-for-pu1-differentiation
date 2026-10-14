function show_results_real_vs_pred_v2(pred_pu1_SS,pred_pu1_prod,figurename)

% pred_pu1_SS is a 1by3 array (from left to right, B cell, pro low, pro
% high)
% pred_pu1_prod is a 1by 3 array

SS_color = '#E03B3B';
prod_color = [0.00 0.75 0.65];
fs = 25;
lw = 3;
ms = 10;
pred_marker_type = 'diamond';
ss_ylim = 1.2;
prod_ylim = 2.5;
load(fullfile(fileparts(pwd), 'Error propagation and Normalization', 'error_prop.mat'));

%%%% categrocial x
cats = categorical({'B cell','pro low','pro high'},'Ordinal',true);
% categories(cats) % debug
desire_order = {'B cell','pro low','pro high'};
desire_cats = reordercats(cats,desire_order);
% categories(desire_cats) % debug
X = double(desire_cats);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% all ss
close all
f = figure('Position',[100 100 800 600]);
axis tight

%%%% real pu1 ss
for ct = 1:3
    curr_conc_range = conc_ratio_all(ct,:);
    errorbar(X(ct),curr_conc_range(2), ...% x, y
        curr_conc_range(2)-curr_conc_range(1),curr_conc_range(3)-curr_conc_range(2), ... % yneg, ypos
        'LineWidth',lw,'Color',SS_color,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',SS_color,'MarkerSize',ms);
    hold on
end

%%%% predicted pu1 ss
plot(X,pred_pu1_SS,'LineStyle','none','linewidth',lw,'Marker',pred_marker_type,'MarkerSize',ms,'MarkerFaceColor','none', 'MarkerEdgeColor',SS_color)
hold on

legend('Data-constrained PU.1-GFP','','','Model-predicted \it{[pu]}', ...
    'Location','northwest')

%%%% save file name
ylim([0 ss_ylim])
xticks(X)
xticklabels(desire_order)
ylabel('Norm. conc.')
set(gca,'FontSize',fs)
exportgraphics(f,strcat(figurename,'pu1_ss','.tif'),'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bcel ss log
f = figure('Position',[100 100 800 600]);
axis tight

%%%% real pu1 ss
for ct = 1:1
    curr_conc_range = conc_ratio_all(ct,:);
    errorbar(X(ct),curr_conc_range(2), ...% x, y
        curr_conc_range(2)-curr_conc_range(1),curr_conc_range(3)-curr_conc_range(2), ... % yneg, ypos
        'LineWidth',lw,'Color',SS_color,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',SS_color,'MarkerSize',ms);
    hold on
end

%%%% predicted pu1 ss
plot(X(1),pred_pu1_SS(1),'LineStyle','none','linewidth',lw,'Marker',pred_marker_type,'MarkerSize',ms,'MarkerFaceColor','none', 'MarkerEdgeColor',SS_color)
hold on

% legend('Model-predicted \it{[pu]}','Data-constrained PU.1-GFP', ...
    % 'Location','northwest')

%%%% save file name
ylim([10^(-4) 10^(-0.5)])
xticks(X(1))
xticklabels(desire_order{1})
ylabel('Norm. conc.')
set(gca,'FontSize',fs,'Yscale','log')
exportgraphics(f,strcat(figurename,'pu1_ss_bcellLog','.tif'),'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% all pu1 prod
f = figure('Position',[100 100 800 600]);
axis tight

%%%% real pu1 prod
for ct = 1:3
    curr_prod_rate_range = prodRate_ratio_all(ct,:);
    errorbar(X(ct),curr_prod_rate_range(2), ...% x, y
        curr_prod_rate_range(2)-curr_prod_rate_range(1),curr_prod_rate_range(3)-curr_prod_rate_range(2), ... % yneg, ypos
        'LineWidth',lw,'Color',prod_color,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',prod_color,'MarkerSize',ms);
    hold on
end

%%%% predicted pu1 prod
plot(X,pred_pu1_prod,'LineStyle','none','linewidth',lw,'Marker',pred_marker_type,'MarkerSize',ms,'MarkerFaceColor','none', 'MarkerEdgeColor',prod_color)
hold on

legend('Data-constrained PU.1 synthesis rate','','','Model-predicted \it[pu]\rm production rate', ...
    'Location','northwest')

%%%% save file name
ylim([0 prod_ylim])
xticks(X)
xticklabels(desire_order)
ylabel('Norm. rate')
set(gca,'FontSize',fs)
exportgraphics(f,strcat(figurename,'pu1_prod','.tif'),'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bcel pu1 ss
f = figure('Position',[100 100 800 600]);
axis tight

%%%% real pu1 prod
for ct = 1:1
    curr_prod_rate_range = prodRate_ratio_all(ct,:);
    errorbar(X(ct),curr_prod_rate_range(2), ...% x, y
        curr_prod_rate_range(2)-curr_prod_rate_range(1),curr_prod_rate_range(3)-curr_prod_rate_range(2), ... % yneg, ypos
        'LineWidth',lw,'Color',prod_color,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',prod_color,'MarkerSize',ms);
    hold on
end

%%%% predicted pu1 prod
plot(X(1),pred_pu1_prod(1),'LineStyle','none','linewidth',lw,'Marker',pred_marker_type,'MarkerSize',ms,'MarkerFaceColor','none', 'MarkerEdgeColor',prod_color)
hold on

% legend('Model-predicted \it[pu]\rm production rate','Data-constrained PU.1 synthesis rate', ...
    % 'Location','northwest')

%%%% save file name
ylim([10^(-2) 10^(-0.5)])
xticks(X(1))
xticklabels(desire_order{1})
ylabel('Norm. rate')
set(gca,'FontSize',fs,'Yscale','log')
exportgraphics(f,strcat(figurename,'pu1_prod_bcellLog','.tif'),'Resolution',300);

