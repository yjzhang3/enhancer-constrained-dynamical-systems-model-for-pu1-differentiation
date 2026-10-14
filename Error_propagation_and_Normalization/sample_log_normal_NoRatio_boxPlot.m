function sample_log_normal_NoRatio_boxPlot(miu_nonmac,sig_nonmac,miu_mac,sig_mac,cellName,ylabelName)

N = 1e6;
nonmac_samples = lognrnd(miu_nonmac,sig_nonmac,[N,1]);
mac_samples = lognrnd(miu_mac,sig_mac,[N,1]);

%%%% plot non-macropahge boxplot
f = figure('Position',[100 100 800 600]);
axis tight
boxplot(nonmac_samples)
set(gca,'FontSize',19)
xlabel(cellName)
ylabel(ylabelName)
if strcmp(ylabelName , 'PU.1-GFP (a.u.)')
    if strcmp(cellName , 'B cell') || strcmp(cellName , 'pro low')
        ylim([0 4.2e4])
    else
        ylim([0 13e4])
    end
elseif strcmp(ylabelName , 'PU.1 synthesis rate (a.u. hr^{-1})')
    if strcmp(cellName , 'B cell') || strcmp(cellName , 'pro low')
        ylim([0 1.7e3])
    else
        ylim([0 4.2e3])
    end
elseif strcmp(ylabelName , 'Cell cycle length (hr)')
    ylim([0 65])
end
exportgraphics(f,strcat('resampled','_',cellName,'_',ylabelName,'.tif'),'Resolution',300);


%%%% plot macropahge boxplot
cellName = 'macrophage';
f = figure('Position',[100 100 800 600]);
axis tight
boxplot(mac_samples)
set(gca,'FontSize',19)
xlabel('macrophage')
ylabel(ylabelName)
if strcmp(ylabelName , 'PU.1-GFP (a.u.)')
    if strcmp(cellName , 'B cell') || strcmp(cellName , 'pro low')
        ylim([0 4.2e4])
    else
        ylim([0 13e4])
    end
elseif strcmp(ylabelName , 'PU.1 synthesis rate (a.u. hr^{-1})')
    if strcmp(cellName , 'B cell') || strcmp(cellName , 'pro low')
        ylim([0 1.7e3])
    else
        ylim([0 6.75e3])
    end
elseif strcmp(ylabelName , 'Cell cycle length (hr)')
    ylim([0 65])
end
exportgraphics(f,strcat('resampled','_',cellName,'_',ylabelName,'.tif'),'Resolution',300);