function sample_log_normal_NoRatio_95CI(miu_nonmac,sig_nonmac,miu_mac,sig_mac,cellName,ylabelName)

N = 1e6;
nonmac_samples = lognrnd(miu_nonmac,sig_nonmac,[N,1]);
mac_samples = lognrnd(miu_mac,sig_mac,[N,1]);

%%%% plot non-macropahge bar
CI = quantile(nonmac_samples, [0.025, 0.975]);
f = figure('Position',[100 100 800 600]);
axis tight
bar(1,mean(nonmac_samples),'FaceColor','none','EdgeColor','k')
hold on
errorbar(1,mean(nonmac_samples),mean(nonmac_samples)-CI(1),CI(2)-mean(nonmac_samples),'LineWidth',3,'Color','k')
set(gca,'FontSize',19)
xticklabels(cellName)
ylabel(ylabelName)
ylim([0 65])
exportgraphics(f,strcat('resampled','_',cellName,'_',ylabelName,'.tif'),'Resolution',300);


%%%% plot macropahge boxplot
cellName = 'macrophage';
f = figure('Position',[100 100 800 600]);
axis tight
CI = quantile(mac_samples, [0.025, 0.975]);
f = figure('Position',[100 100 800 600]);
axis tight
bar(1,mean(mac_samples),'FaceColor','none','EdgeColor','k')
hold on
errorbar(1,mean(mac_samples),mean(mac_samples)-CI(1),CI(2)-mean(mac_samples),'LineWidth',3,'Color','k')
set(gca,'FontSize',19)
xticklabels(cellName)
ylabel(ylabelName)
ylim([0 65])
exportgraphics(f,strcat('resampled','_',cellName,'_',ylabelName,'.tif'),'Resolution',300);