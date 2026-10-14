function plotonepar_wptName(x,s,rangePar,f,ivar,PlotOptions,coli,axNum)
% ivar is the ith variable to be plotted (index of the species, GATA2 is
% the second species therefores 2)

if rangePar ~= 0
    % 1-Par Bifurcation diagram 
    ns = length(s);
    for j=1:ns-1 % this is basically plotting all "segments" of lines between every two found points 
    %     if j > 1
    %         hvis = 'off';
    %     end
        A = f(:,s(j).index:s(j+1).index); % extract eigenvalue matrix for that segment
        mx = max(real(A), [], 1);
        allNeg = all(mx < 0);
        if allNeg == 1 % check stability by looking at the *second point* (not segment??) on the currently plotted segment
             plot(axNum,rangePar(s(j).index:s(j+1).index), x(ivar,s(j).index:s(j+1).index),...
            'Color',coli,'LineStyle','-','LineWidth',PlotOptions.lw)
    %     'DisplayName',lname,'HandleVisibility',hvis); % plot Gata2 species over parameter 1 (BMP)
        else 
            plot(axNum,rangePar(s(j).index:s(j+1).index), x(ivar,s(j).index:s(j+1).index),...
            'Color',coli,'LineStyle','-','LineWidth',PlotOptions.lw)
            hold on

            % extract violating columns (where max real part >= 0)
            badCols = A(:, mx >= 0);
            badCols_idx = find(mx >= 0);
            plot(axNum,rangePar(badCols_idx), x(ivar,badCols_idx),...
            'Color','y','LineStyle',':','LineWidth',PlotOptions.lw)
        end
        hold on;
        plot(axNum,rangePar(s(j).index),x(ivar,s(j).index),'Marker','*','MarkerSize',10,'MarkerEdgeColor','r') % mark the singulariy point
    %     ,'HandleVisibility','off'); % plot the LP or HP
        hold on
        text(axNum,rangePar(s(j).index),x(ivar,s(j).index),s(j).label,'FontSize',12,'HorizontalAlignment','center') % label the point type
        hold on
    end
    hold on; axis square; box on;
    % xlim([0 max(rangePar)]);
    % ylim([0 max(x(ivar,:))])
    set(gca,'XScale',PlotOptions.xsc,'YScale',PlotOptions.ysc,...
            'XTickLabelRotation', 0,...
            'TickLength',PlotOptions.tckL,'FontSize',PlotOptions.Fs,...
            'FontName','Arial');
    %         'XTick',PlotOptions.xtk,'YTick',PlotOptions.ytk,...,
else
    return
end