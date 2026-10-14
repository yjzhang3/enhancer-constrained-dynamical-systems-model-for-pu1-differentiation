% why not plot the kueh paper proposed phenemological growth rate, and then
% plot the bar graph of grwoth rate (converted from cell cycle), and point
% out the contradiction?
clear
clc
close all
% close all
fs = 25; % overall fontsize
fs_small = 15; % legend fontsize
lw = 3; % linewidth
ms = 15; % markersize
mk = 'o'; % marker

%% plot kueh model "cell division rate" model
pu1_totest_range = 0:0.01:1.2;
celldiv_rate_func = @(Cmin,Cmax,Kc,N,x) Cmin + (Cmax-Cmin)*Kc^N/(Kc^N+x^N);
Cmin = 0.02;
Cmax = 0.1;
Kc = 20000;
N = 20;
% D = 0.02;
p_ss = 37500;

dilution_rate_range = zeros(size(pu1_totest_range));
for n = 1:length(pu1_totest_range)
    dilution_rate_range(n) = celldiv_rate_func(Cmin,Cmax,Kc,N,pu1_totest_range(n)*p_ss)/(celldiv_rate_func(Cmin,Cmax,Kc,N,p_ss));
end

%% plot pu1-dependent growth rate (as an inset) 
f = figure('Position',[100 100 800 600]);
axis tight

plot(pu1_totest_range,dilution_rate_range,'LineWidth',lw,'Color','k','LineStyle','--')
xlabel('\it{[pu]}')
ylabel('Normalized cell division rate')
ylim([0.5 5])
set(gca,'FontSize',fs)
exportgraphics(f,'kueh cell division normalized rate.tif','Resolution',300);


%% plot experimental cell cycle normalized by macrophage removal rate
f = figure('Position',[100 100 800 600]);
axis tight
load(fullfile(fileparts(pwd), 'Error propagation and Normalization', 'error_prop.mat'));
pu1_bcell_conc_range = conc_ratio_all(1,:); % skip the median of ratio
pu1_ProLow_conc_range = conc_ratio_all(2,:);
pu1_ProHigh_conc_range = conc_ratio_all(3,:);
pu1_conc_all = [pu1_bcell_conc_range;pu1_ProLow_conc_range;pu1_ProHigh_conc_range];

cellcycle_real_Bcell = cellcycle_length_ratio_all(1,:);
cellcycle_real_ProLow = cellcycle_length_ratio_all(2,:);
cellcycle_real_ProHigh = cellcycle_length_ratio_all(3,:);
cellcycle_real_all = [cellcycle_real_Bcell;cellcycle_real_ProLow;cellcycle_real_ProHigh];

% a dummy point
errorbar(0.8,0.3,.02,.02,.02,.02,'LineStyle','none','Marker',mk,'MarkerFaceColor','k','MarkerEdgeColor','none','Color','k','MarkerSize',ms,'LineWidth',lw)
% h.Visible = 'off';

hold on

% name_arr = {'C_{B cell}','C_{pro low}','C_{pro high}','C_{macrophage}'};
name_arr = {'B cell','pro low','pro high','macrophage'};
color_arr = {'#FF5800','#BFBFBF','#808080','#2626B2'};

text_vertical_delta = 0.2;
for y = 1:3
    arr_cc = cellcycle_real_all(y,:);
    arr_conc = pu1_conc_all(y,:);
    h = errorbar(arr_conc(2),arr_cc(2), ... % x , y
        arr_cc(1)-arr_cc(2),arr_cc(2)-arr_cc(3), ... % yneg, ypos
        arr_conc(2)-arr_conc(1),arr_conc(3)-arr_conc(2), ... % xneg,xpos
        'LineStyle','none','Marker',mk,'MarkerFaceColor',color_arr{y},'MarkerEdgeColor','none','Color',color_arr{y},'MarkerSize',ms,'LineWidth',lw);
    text(arr_conc(2),arr_cc(2)+text_vertical_delta,name_arr{y},'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fs_small)
    hold on
end
y = 4;
plot(1,1,'LineStyle','none','Marker',mk,'MarkerFaceColor',color_arr{y},'MarkerEdgeColor','none','Color',color_arr{y},'MarkerSize',ms,'LineWidth',lw);
text(1,1+text_vertical_delta,name_arr{y},'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fs_small)


%% plot cell cycle length
cellcycle_range = log(2)./dilution_rate_range;
plot(pu1_totest_range,cellcycle_range,'LineWidth',lw,'Color','k','LineStyle','-')
hold on

%% overall
legend(['Experimental cell cycle length',newline, ''],'','','','',['\it{[pu]}\rm-predicted cell cycle length',newline, ''],'Location','Northwest','FontSize',fs_small)
xlabel('\it{[pu]}')
ylabel('Norm. cell cycle length')
xlim([0 1.2])
ylim([0 1.5])
ax = gca;
ax.Box = 'on';
set(gca,'FontSize',fs)
exportgraphics(f,'contradicting pu1-dependent cc v2.tif','Resolution',300);