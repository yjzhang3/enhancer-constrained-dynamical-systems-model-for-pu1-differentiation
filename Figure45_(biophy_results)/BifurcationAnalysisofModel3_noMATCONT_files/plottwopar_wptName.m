function plottwopar_wptName(x_ap2f,s_ap2f,f_ap2f,x_ap2b,s_ap2b,f_ap2b,updated_apn,PlotOptions)

%% test stability region
stabf = test_stability(f_ap2f(:,2:end-1));
stabb = test_stability(f_ap2b(:,2:end-1));

%% extract values to be plotted on horizontal and vertical axis
num_x = length(x_ap2f(:,1))-2; % number of state variables
range_ap12f_hor = x_ap2f(num_x+2,:);
range_ap12f_vert = x_ap2f(num_x+1,:);
range_ap12b_hor = x_ap2b(num_x+2,:);
range_ap12b_vert = x_ap2b(num_x+1,:);

%% plot forward with point name
ns = length(s_ap2f); % number of line segment
for j = 1:ns-1
    plot(range_ap12f_hor(s_ap2f(j).index:s_ap2f(j+1).index), range_ap12f_vert(s_ap2f(j).index:s_ap2f(j+1).index),...
        'Color',PlotOptions.colorf,'LineWidth',PlotOptions.lw)
    hold on
    plot(range_ap12f_hor(s_ap2f(j).index),range_ap12f_vert(s_ap2f(j).index),'Marker','*','MarkerSize',10,'MarkerEdgeColor','r') % mark the singulariy point
    hold on
    text(range_ap12f_hor(s_ap2f(j).index),range_ap12f_vert(s_ap2f(j).index),s_ap2f(j).label,'FontSize',12) % label the point type
end
hold on
%% plot backward with point name
ns = length(s_ap2f); % number of line segment
clear j
for j = 1:ns-1
    plot(range_ap12b_hor(s_ap2f(j).index:s_ap2b(j+1).index), range_ap12b_vert(s_ap2f(j).index:s_ap2b(j+1).index),...
        'Color',PlotOptions.colorb,'LineWidth',PlotOptions.lw)
    hold on
    plot(range_ap12b_hor(s_ap2f(j).index),range_ap12b_vert(s_ap2f(j).index),'Marker','*','MarkerSize',10,'MarkerEdgeColor','r') % mark the singulariy point
    hold on
    text(range_ap12b_hor(s_ap2f(j).index),range_ap12b_vert(s_ap2f(j).index),s_ap2b(j).label,'FontSize',12) % label the point type
end

hold on; axis square; box on;
xlim(PlotOptions.xL);
ylim(PlotOptions.yL);
set(gca,'XScale',PlotOptions.xsc,'YScale',PlotOptions.ysc,...
        'XTickLabelRotation', 0,...
        'TickLength',PlotOptions.tckL,'FontSize',PlotOptions.Fs,...
        'FontName','Arial');
%         'XTick',PlotOptions.xtk,'YTick',PlotOptions.ytk,...,
xlabel(updated_apn{2})
ylabel(updated_apn{1})
title(strcat(updated_apn{1},{' '},'vs',{' '},updated_apn{2},{' '},'two-par'))