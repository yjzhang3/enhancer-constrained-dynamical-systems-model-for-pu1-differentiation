function [good_pars,D7_prob_bcell,minerr,best_pars] = predict_cebpa_PU1CEfeedback(objFuncHandle,selec,final_fn,nuc_destab_location)
% Subdirectory name
subdir_fig = 'outputs_figs';
subdir_mat = 'outputs_mats';

% Ensure subdirectory exists (creates if not)
if ~exist(subdir_fig, 'dir')
    mkdir(subdir_fig);
end

if ~exist(subdir_mat, 'dir')
    mkdir(subdir_mat);
end

%% read file 
if ~isfile(final_fn)
    disp('file does not exist')
    good_pars = 0;
    D7_prob_bcell = 0;
    minerr = NaN;
    return
else
    n_datap = length(selec); % number of data points to fit
    A = readmatrix(final_fn);
end

% find degenerate sets of pars just by looking at error
arri = A(:,end);
sorted_arri = sort(arri);
minerr = sorted_arri(1); % pick the smallest or second smallest?
minind = find(arri==minerr);
% [minerr,minind] = min(arr_to_sort);
tol_err = 0.005;
good_idx = find(A(:,end)<= minerr*(1+tol_err));

% find degenerate sets of pars just by looking at euclidean difference
% between all pars
min_params = A(minind,1:end-1);
tol_params = 0.05;
distances = vecnorm(A(good_idx,1:end-1) - min_params, 2, 2); % Euclidean
better_idx = good_idx(distances < tol_params);

%% other logistics
nvars = length(A(1,:))-1;
changePU1conc = 0;
notxt = erase(final_fn,".txt"); 
% notxt = erase(notxt,'merged_files_0924/');
% orig_name = strrep(notxt,"_"," ");
% objFuncHandle = @objective_with_cebpa_runx1_general_alt1_PU1CEfeedback;

%% extract all good fits
good_pred = zeros(length(better_idx),n_datap);
good_pars = zeros(length(better_idx),nvars);
for yt = 1:length(better_idx)
    [err,sim_d_final,real_d_final] = objFuncHandle(selec,A(better_idx(yt),1:end-1),changePU1conc,0,0,nuc_destab_location);
    disp(err)
    good_pred(yt,:) = sim_d_final;
    good_pars(yt,:) = A(better_idx(yt),1:end-1);
end
row_num = 1;
best_pred = good_pred(row_num,:);
best_pars = good_pars(row_num,:);

%% plot prediction vs. experimental

load('Leddin_real_d_std.mat','real_d_std')
real_d_std = [real_d_std,10];
all_data_label = {'D2/D1','D3/D1','D4/D1','D5/D1','D6/D1','D7/D1','F2/E1','F3/E1','F4/E1','E2/E1','E3/E1'};

f=figure();
% errorbar(1:n_datap,mean(best_pred),std(best_pred),'LineStyle','none','Marker','d','MarkerSize',20,'Color','m')
plot(1:n_datap,best_pred,'LineStyle','none','Marker','d','MarkerSize',20,'Color',[0.00 0.75 0.65])
hold on
errorbar(1:n_datap,real_d_final,real_d_std(selec)./2556.3,'LineStyle','none','Marker','o','MarkerSize',20,'Color','k')
set(gca, 'YScale', 'linear')
xticks([1:n_datap])
xticklabels(all_data_label(selec))
ylabel('Ratio in gene expression')
set(gca,'FontSize',15)
destab_id = num2str(nuc_destab_location);
final_id = strrep(destab_id," ","");
% legend(strcat('Model','(',final_id,')'),'Experiment','Location','Northwest')
legend('Model','Experiment','Location','Northwest')
% title(strrep(orig_name, '_', ' '))
title(final_id)
saveas(f,fullfile(subdir_fig,strcat('prediction_',notxt,'.jpg')))

%%
% figure('position',[344,41,1178,994])
% f=figure('Units', 'normalized', 'Position', [0.05, 0.05, 0.8, 1.2]);
% ax1 = axes('Position',[0.2, 0.4, 0.75, 0.55]); %[left, bottom, width, height]
% % subplot(2,1,1)
% colororder([0.5,0.5,0.5])
% 
% % energy_table = [G_NUC,G_PU1_14,G_PU1_12,G_PU1_10,G_CEBP,G_RUNX,G_RNAP, ... %1-7
% %     G_NUC_PU1,G_NUC_CEBP,G_NUC_RUNX1,... %8-10 
% %     G_PU1_RNAP,... % 11
% %     G_RUNX_PU1,G_RUNX_CEBP,G_CEBP_PU1]; ... % 12-14 
% 
% allnames_G = {'G_{nuc}','G_{PU114kb}','G_{PU112kb}','G_{PU110kb}','G_{CEBP}','G_{RNAP}',...
%     'G_{PU1-nuc}','G_{CEBP-nuc}','G_{RNAP-PU1}','G_{RNAP-CEBP}',...
%     'G_{PU1-C/EBPa}'};
% X = categorical(allnames_G);
% X = reordercats(X,allnames_G);
% Y = best_pars(1:length(allnames_G));
% barh(ax1,X,Y);
% set(ax1,'FontSize',20,'Xscale','linear')
% 
% hold on
% ax2 = axes('Position',[0.2, 0.1, 0.75, 0.25]); %[left, bottom, width, height]
% allnames_other = {'rescaling factor','[pu1]_{B cell}','[pu1]_{Pro Low}','[pu1]_{Pro High}',...
%     'Amin','Amax','N','K'};
% X = categorical(allnames_other);
% X = reordercats(X,allnames_other);
% Y = best_pars(length(allnames_G)+1:end);
% bar(ax2,X,Y);
% set(ax2,'FontSize',20,'YMinorTick','on','Yscale','linear')
% ylim([0 6])
% % title(all_data_label(selec))
% % xline(1,'LineStyle','--','Color','r')
% % sgtitle(strrep(orig_name, '_', ' '))
% sgtitle(final_id)
% 
% % set(gca,'FontSize',30,'XMinorTick','on','Xscale','linear')
% % title(strrep(orig_name, '_', ' '))
% yline(1,'LineStyle','--','Color','r')
% saveas(f,fullfile('outputs_figs',strcat('reportG',notxt,'.jpg')))

%% 
save(fullfile(subdir_mat,strcat(notxt,'.mat')),'best_pars')