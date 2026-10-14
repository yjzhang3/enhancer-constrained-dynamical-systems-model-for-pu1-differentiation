clear
clc
close all
parsType = {'Amin','K'};
% for ct = 1:3
for ct = 3
    % for p = 1:length(parsType)
    for p = 2
        loadFigs(ct,parsType{p})
    end
end



function loadFigs(ct,parsType)
%% housekeeping
% clear
% clc
close all
bcellcolor = '#FF5800';
prolowcolor = '#BFBFBF';
prohighcolor = 'k';
color_scheme = {bcellcolor,prolowcolor,prohighcolor};
fa = 0.6; % facealpha of constraint patch (not growth rate patch)
load('error_prop_0216.mat')
% ct = 3;
% parsType = 'K';
outputdir = strcat('Results',num2str(ct),'/',parsType);
if ~exist(outputdir, 'dir')
    mkdir(outputdir);
end

%% load error plot
close all
fn2 = strcat('CT=',num2str(ct),'_cebphill_{',parsType,'}_all_err_vs_rangealpha.fig');
% ct = str2double(fn2(4));
fnew = openfig(fn2); 
set(fnew, 'Position', [268,465,751,574]);
ylabel('Model error')

ax2 = gca;
set(ax2,'FontSize',16)

% find minimum error
% Grab every plotted object that has YData (lines, scatter, errorbar, etc.)
h = findobj(ax2, '-property', 'YData');

allX = [];
allY = [];
for i = 1:numel(h)
    allX = [allX; h(i).XData(:)];
    allY = [allY; h(i).YData(:)];
end

[minY, idx] = min(allY);
minX = allX(idx);

hold on;
plot(minX, minY, 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);   % optional: mark the point
text(minX, minY, sprintf('  min = %.3g', minY), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
legend('off')

if contains(fn2,'Amin')
    exL = 21;
else
    exL = 18;
end
currentTitle = ax2.Title.String;
ax2.Title.String = currentTitle(1:min(exL, length(currentTitle)));

hold off;
exportgraphics(ax2,fullfile(outputdir,strcat(fn2(1:end-4),'_final','.tif')),'Resolution','300')

%% load bifurcation diagram, overlay with color patch
fn = strcat('[x]_vs_cebphill_{',parsType,'}_ct=',num2str(ct),'.fig');
ct = str2double(fn(end-4));
fnew=openfig(fn); 

axAll = findall(fnew, 'Type', 'axes');

% findobj often returns axes in reverse creation order (last-created first),
% so sort by their horizontal position to reliably get left-to-right order
positions = vertcat(axAll.Position);   % Nx4 matrix, works for N=1 too
[~, idx] = sort(positions(:,1));   % sort by x (left edge) position
axAll = axAll(idx);
ax1 = axAll(1);   % leftmost subplot

hold on
xl = xlim;   % get the current x-axis limits [xmin xmax] from the open figure
y0 = conc_ratio_all(ct,1);      % baseline y-value where the patch starts (adjust to your needs)
h = conc_ratio_all(ct,3) - conc_ratio_all(ct,1);

hold on;
p = patch(ax1,[xl(1) xl(2) xl(2) xl(1)], [y0 y0 y0+h y0+h], [0.8 0.8 0.8], ...
    'FaceAlpha', fa, 'FaceColor',color_scheme{ct},'EdgeColor','none');
xlim(xl)
uistack(p, 'bottom');   % optional: push patch behind existing data points/lines

% xline(ax1,minX,'Color','b','LineStyle','--','LineWidth',4)
h = findobj(ax1, '-property', 'YData');

allX2 = [];
allY2 = [];
for i = 1:numel(h)
    allX2 = [allX2; h(i).XData(:)];
    allY2 = [allY2; h(i).YData(:)];
end
xtoplot=allX2(find(allX2==minX));
ytoplot=allY2(find(allX2==minX));

plot(ax1,xtoplot(1),ytoplot(1),'LineStyle','none','Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10)

% yline(ax1,ytoplot(1),'Color','b','LineStyle','--','LineWidth',4)
% text(ax1,xtoplot(1)-0.05,0.05+ytoplot(1),sprintf('par=%.4f,[pu]=%.4f',minX,ytoplot(1)),'FontSize',15)

hold off;
exportgraphics(fnew,fullfile(outputdir,strcat(fn(1:end-4),'_final','.tif')),'Resolution','300')


end