function test_Hill_property_cebpa_general_PU1CEfeedback(objFuncHandle,final_fn,nuc_destab_location)
if ~isfile(final_fn)
    disp('file does not exist')
    return
else
    load(final_fn,'best_pars')
end

destab_id = num2str(nuc_destab_location);
final_id = strrep(destab_id," ","");

% subdir_fig = 'outputs_figs_ss';
% % Ensure subdirectory exists (creates if not)
% if ~exist(subdir_fig, 'dir')
%     mkdir(subdir_fig);
% end

% load('cell_cycle_pars_Kueh.mat','C_max','C_min','Kc','N','D')
% load('cell_cycle_pars_true.mat','C_max','C_min','D') % if want to test biologially cmax, cmin, and D
fn = erase(final_fn,".mat");
fn = erase(fn,'outputs_mats/');
fn = "";
linestyle_name = {'-','--',':','-o','->','-v','-square'};
lw = 3; % linewidth
ms = 12; % marker size
fs = 19; % fontsize
fa = 0.6; % facealpha of constraint patch (not growth rate patch)
% color for constraint patch
bcellcolor = '#FF5800';
prolowcolor = '#BFBFBF';
prohighcolor = 'k';
maccolor = '#2626B2';

% color and style for removal rate linear line
removalcolor = [0.70 0.40 0.85];
prodcolor = [0.00 0.75 0.65];
removalls = '-';

% color of intersection
color_mk_intersect = '#E03B3B';

% upper limit of plot
xlim_conc = 1.2;
ylim_prod_rate = 2.5;


% objFuncHandle = @objective_with_cebpa_runx1_general_alt1_PU1CEfeedback;
RNAPFuncHandle = @RNAP_occ_prob_cebpa_general;

%% just test relaly quick to see if the best pars is indeed correct
p = best_pars(1,:);
[final_er,sim_d_final,real_d_final,energy_table,concentration_table,concentration_table_Bcell,...
    concentration_table_ProLow,concentration_table_ProHigh,all_sublist_str,...
    norm_D7_prob_bcell,norm_D7_prob_ProLow,norm_D7_prob_ProHigh,p_cebp_hill,pu1_prod_slope_celltype] = objFuncHandle([1:11],p,0,0,0,nuc_destab_location);
disp(final_er)
D7_valid_sublist = all_sublist_str.('D7');
model_predicted_Bcell_pu1 = concentration_table_Bcell(2);
model_predicted_Bcell_cebpa = concentration_table_Bcell(3);

model_predicted_ProLow_pu1 = concentration_table_ProLow(2);
model_predicted_ProHigh_pu1 = concentration_table_ProHigh(2);
model_predicted_ProLow_cebpa = concentration_table_ProLow(3);
model_predicted_ProHigh_cebpa = concentration_table_ProHigh(3);

disp([model_predicted_Bcell_pu1 model_predicted_ProLow_pu1 model_predicted_ProHigh_pu1])

% save(sprintf('cebpa_PU1CEfeedback_prediction_PU1SS_N%.0f.mat',p_cebp_hill(3)), ...
%     'norm_D7_prob_bcell','norm_D7_prob_ProLow','norm_D7_prob_ProHigh',...
%     'model_predicted_Bcell_pu1','model_predicted_ProLow_pu1','model_predicted_ProHigh_pu1','pu1_prod_slope_celltype',...
%     'model_predicted_Bcell_cebpa','model_predicted_ProLow_cebpa','model_predicted_ProHigh_cebpa')

%% ss production rate
D7_ss_prod_rate = sim_d_final(6);

%% Range of [PU1] concentration and synthesis rate in which we find intersection (color patch)
load(fullfile(parentdir(pwd, 1), 'Error propagation and Normalization', 'error_prop.mat'));
pu1_bcell_conc_range = conc_ratio_all(1,[1,3]); % skip the median of ratio because I'm plotting patch
pu1_ProLow_conc_range = conc_ratio_all(2,[1,3]);
pu1_ProHigh_conc_range = conc_ratio_all(3,[1,3]);

pu1_bcell_prod_rate_range = prodRate_ratio_all(1,[1,3]);
pu1_ProLow_prod_rate_range = prodRate_ratio_all(2,[1,3]);
pu1_ProHigh_prod_rate_range = prodRate_ratio_all(3,[1,3]);

f1 = [1 2 3 4];

% what if we combine production and [pu1] range?
color_patch_pu1_bcell_overall = [ pu1_bcell_conc_range(1) pu1_bcell_prod_rate_range(1); pu1_bcell_conc_range(2) pu1_bcell_prod_rate_range(1); pu1_bcell_conc_range(2) pu1_bcell_prod_rate_range(2); pu1_bcell_conc_range(1) pu1_bcell_prod_rate_range(2)]; % vertices of color patch
color_patch_pu1_ProLow_overall = [ pu1_ProLow_conc_range(1) pu1_ProLow_prod_rate_range(1); pu1_ProLow_conc_range(2) pu1_ProLow_prod_rate_range(1); pu1_ProLow_conc_range(2) pu1_ProLow_prod_rate_range(2); pu1_ProLow_conc_range(1) pu1_ProLow_prod_rate_range(2)]; % vertices of color patch
color_patch_pu1_ProHigh_overall = [ pu1_ProHigh_conc_range(1) pu1_ProHigh_prod_rate_range(1); pu1_ProHigh_conc_range(2) pu1_ProHigh_prod_rate_range(1); pu1_ProHigh_conc_range(2) pu1_ProHigh_prod_rate_range(2); pu1_ProHigh_conc_range(1) pu1_ProHigh_prod_rate_range(2)]; % vertices of color patch

%% range of growth rate
rem_rate_constant_Bcell = remRate_constant_ratio_all(1,2);
rem_rate_constant_Bcell_Q1 = remRate_constant_ratio_all(1,1);
rem_rate_constant_Bcell_Q3 = remRate_constant_ratio_all(1,3);

rem_rate_constant_ProLow = remRate_constant_ratio_all(2,2);
rem_rate_constant_ProLow_Q1 = remRate_constant_ratio_all(2,1);
rem_rate_constant_ProLow_Q3 = remRate_constant_ratio_all(2,3);

rem_rate_constant_ProHigh = remRate_constant_ratio_all(3,2);
rem_rate_constant_ProHigh_Q1 = remRate_constant_ratio_all(3,1);
rem_rate_constant_ProHigh_Q3 = remRate_constant_ratio_all(3,3);

rem_rate_constant_Mac = 1;
% 
% model_predicted_Bcell_rem_rate_range = model_predicted_Bcell_pu1*remRate_constant_ratio_all(1,:);
% model_predicted_ProLow_rem_rate_range = model_predicted_ProLow_pu1*remRate_constant_ratio_all(2,:);
% model_predicted_ProHigh_rem_rate_range = model_predicted_ProHigh_pu1*remRate_constant_ratio_all(3,:);
% 
% save(sprintf('cebpa_PU1CEfeedback_prediction_PU1SS_N%.0f.mat',p_cebp_hill(3)),...
%     'model_predicted_Bcell_rem_rate_range','model_predicted_ProLow_rem_rate_range','model_predicted_ProHigh_rem_rate_range','-append')

%% relative expression of cebpa and runx1 in b cell
% cebp_cell_type_specific_conc_bcell = [4/46*1.42/4.77,4/43,model_predicted_Bcell_cebpa];
% runx1_cell_type_specific_conc_bcell = [model_predicted_Bcell_runx1,1/100];
cebp_cell_type_specific_conc_bcell = [model_predicted_Bcell_cebpa];
% runx1_cell_type_specific_conc_bcell = [];

% all possible combination
% [m, n] = meshgrid(cebp_cell_type_specific_conc_bcell, runx1_cell_type_specific_conc_bcell');
% bcell_all_possbility = [reshape(m, [], 1), reshape(n, [], 1)]; 
bcell_all_possbility = cebp_cell_type_specific_conc_bcell;
bcell_all_possibility_legend = cell(length(bcell_all_possbility(:,1)),1);
bcell_all_possibility_linestyle = cell(length(bcell_all_possbility(:,1)),1);
for i = 1:length(bcell_all_possbility(:,1))
    bcell_all_possibility_legend{i} = sprintf('B cell, \\it[ce]\\rm = %.3f',bcell_all_possbility(i,1));
    bcell_all_possibility_linestyle{i} = linestyle_name{i};
end

%% relative expression of cebpa and runx1 in prog low
% cebp_cell_type_specific_conc_pro = [1.42/4.77,0.038/0.37,model_predicted_Bcell_cebpa*46/4];
% runx1_cell_type_specific_conc_pro = [0.18/0.22];
cebp_cell_type_specific_conc_prolow = [model_predicted_ProLow_cebpa];
% runx1_cell_type_specific_conc_prolow = [model_predicted_Pro_runx1];

% all possible combination
% [m, n] = meshgrid(cebp_cell_type_specific_conc_prolow, runx1_cell_type_specific_conc_prolow');
% prolow_all_possibility = [reshape(m, [], 1), reshape(n, [], 1)]; 
prolow_all_possibility = cebp_cell_type_specific_conc_prolow;
prolow_all_possibility_legend = cell(length(prolow_all_possibility(:,1)),1);
prolow_all_possibility_linestyle = cell(length(prolow_all_possibility(:,1)),1);

for i = 1:length(prolow_all_possibility(:,1))
    prolow_all_possibility_legend{i} = sprintf('pro low, \\it[ce]\\rm = %.3f',prolow_all_possibility(i,1));
    prolow_all_possibility_linestyle{i} = linestyle_name{i};
end

% color scheme of plot
% load('best_fit_model_prediction_colors.mat','hexColor')
% rgbColor = hex2rgb(hexColor);
% rgbColor = [rgbColor(2:3,:);[1,0,1]];

%% relative expression of cebpa and runx1 in prog high
% cebp_cell_type_specific_conc_pro = [1.42/4.77,0.038/0.37,model_predicted_Bcell_cebpa*46/4];
% runx1_cell_type_specific_conc_pro = [0.18/0.22];
cebp_cell_type_specific_conc_prohigh = [model_predicted_ProHigh_cebpa];
% runx1_cell_type_specific_conc_prohigh = [model_predicted_Pro_runx1];

% all possible combination
% [m, n] = meshgrid(cebp_cell_type_specific_conc_prohigh, runx1_cell_type_specific_conc_prohigh');
% prohigh_all_possibility = [reshape(m, [], 1), reshape(n, [], 1)]; 
prohigh_all_possibility = cebp_cell_type_specific_conc_prohigh;
prohigh_all_possibility_legend = cell(length(prohigh_all_possibility(:,1)),1);
prohigh_all_possibility_linestyle = cell(length(prohigh_all_possibility(:,1)),1);

for i = 1:length(prohigh_all_possibility(:,1))
    prohigh_all_possibility_legend{i} = sprintf('pro high, \\it[ce]\\rm = %.3f',prohigh_all_possibility(i,1));
    prohigh_all_possibility_linestyle{i} = linestyle_name{i};
end

% color scheme of plot
% load('best_fit_model_prediction_colors.mat','hexColor')
% rgbColor = hex2rgb(hexColor);
% rgbColor = [rgbColor(2:3,:);[1,0,1]];

%% all possible pu1 level
pu1_totest_range = [0:0.001:xlim_conc];

%% bcell
f = figure('Position',[100 100 800 600]);
axis tight
plot_option = 'log';

% %%%%%%% plot overaall color patch
patch('Faces',f1,'Vertices',color_patch_pu1_bcell_overall,'FaceColor',bcellcolor,'EdgeColor','none','FaceAlpha',fa*0.75)
hold on
lgd = "Data-constrained B cell";
hold on

% %%%%%%% plot production rate
TR_range = zeros(size(pu1_totest_range));

for m = 1:length(bcell_all_possbility(:,1))
    cebpa = bcell_all_possbility(m,1);
    % runx1 = bcell_all_possbility(m,2);
    for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
        
        D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        
        D7_curr_prod_rate = D7_prob;

        TR_range(i) = D7_curr_prod_rate/D7_ss_prod_rate;
    end
    
    plot(pu1_totest_range,TR_range,bcell_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
    hold on

    cell_name = bcell_all_possibility_legend{m};
    lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
end

% %%%%%%% generate growth curve
% tau_Bcell = 25.5;
% tau_Bcell_Q1 = 17.7;
% tau_Bcell_Q3 = 41.2;
min_removal_rate_slope = rem_rate_constant_Bcell_Q1;
max_removal_rate_slope = rem_rate_constant_Bcell_Q3;

removal_range = rem_rate_constant_Bcell*pu1_totest_range;
min_removal_rate = min_removal_rate_slope*pu1_totest_range;
max_removal_rate = max_removal_rate_slope*pu1_totest_range;

plot(pu1_totest_range,removal_range,'LineWidth',lw,'Color',removalcolor,'LineStyle',removalls)
hold on
plot_range_growthrate_patch(pu1_totest_range,max_removal_rate,min_removal_rate,removalcolor)
hold on

% %%%%%% find actual steady state
[~,ind] = min(abs(TR_range-removal_range));
concise_ind = extract_intersection_index(ind,40);
true_Bcell_pu1_ss = pu1_totest_range(concise_ind);
true_Bcell_pu1_prod = TR_range(concise_ind);

% %%%%%% mark itnersection
% xline(model_predicted_Bcell_pu1,'Color','r','LineWidth',2)
% yline(D7_prob_bcell/D7_ss_prod_rate,'Color','r','LineWidth',2)
% plot(true_Bcell_pu1_ss,true_Bcell_pu1_prod,'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
plot(model_predicted_Bcell_pu1,norm_D7_prob_bcell,'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
% text(model_predicted_Bcell_pu1,norm_D7_prob_bcell,"WT PU.1 production in B cell",'FontSize',12,'Color','r','HorizontalAlignment','left')

% %%%%%% together
lgd = [lgd,"Remv. rate in B cell","",""];
legend(lgd,'location','northwest')
xlabel('\it{[pu]}')
ylabel('Normalized Rate')
set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)

if contains(plot_option,'log')
    plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
    plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
else
    plot_ub_x = [0,pu1_totest_range(end)];
    plot_ub_y = [0,ylim_prod_rate];
end
xlim(plot_ub_x)
ylim(plot_ub_y)

% title(strrep(fn, '_', ' '))
% title(final_id)
% saveas(f,fullfile('outputs_figs',strcat(fn,'_Bcellonly','_',plot_option,'.tif')))
ax = gca;
ax.Box = 'on';
%%exportgraphics(f,fullfile('outputs_figs',strcat(fn,'_Bcellonly','_',plot_option,'.tif')),'Resolution',300);




%% bcell log
% f = figure('Position',[100 100 800 600]);
% axis tight
% plot_option = 'log';
% 
% % %%%%%%% plot overaall color patch
% patch('Faces',f1,'Vertices',color_patch_pu1_bcell_overall,'FaceColor',bcellcolor,'EdgeColor','none','FaceAlpha',fa*0.75)
% hold on
% lgd = "Data-constrained B cell";
% hold on
% 
% % %%%%%%% plot production rate
% TR_range = zeros(size(pu1_totest_range));
% 
% for m = 1:length(bcell_all_possbility(:,1))
%     cebpa = bcell_all_possbility(m,1);
%     % runx1 = bcell_all_possbility(m,2);
%     for i = 1:length(pu1_totest_range)
%         pu1 = pu1_totest_range(i);
%         totest_concentration_table = [1,pu1,cebpa,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
% 
%         D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
% 
%         D7_curr_prod_rate = D7_prob;
% 
%         TR_range(i) = D7_curr_prod_rate/D7_ss_prod_rate;
%     end
% 
%     plot(pu1_totest_range,TR_range,bcell_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
%     hold on
% 
%     cell_name = bcell_all_possibility_legend{m};
%     lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
% end
% 
% % %%%%%%% generate growth curve
% % tau_Bcell = 25.5;
% % tau_Bcell_Q1 = 17.7;
% % tau_Bcell_Q3 = 41.2;
% min_removal_rate_slope = rem_rate_constant_Bcell_Q1;
% max_removal_rate_slope = rem_rate_constant_Bcell_Q3;
% 
% removal_range = rem_rate_constant_Bcell*pu1_totest_range;
% min_removal_rate = min_removal_rate_slope*pu1_totest_range;
% max_removal_rate = max_removal_rate_slope*pu1_totest_range;
% 
% plot(pu1_totest_range,removal_range,'LineWidth',lw,'Color',removalcolor,'LineStyle',removalls)
% hold on
% plot_range_growthrate_patch(pu1_totest_range,max_removal_rate,min_removal_rate,removalcolor)
% hold on
% 
% % %%%%%% find actual steady state
% ind = find(abs(TR_range-removal_range)<0.00001);
% concise_ind = extract_intersection_index(ind,40);
% true_Bcell_pu1_ss = pu1_totest_range(concise_ind);
% true_Bcell_pu1_prod = TR_range(concise_ind);
% 
% % %%%%%% mark production and growth rate
% % xline(model_predicted_Bcell_pu1,'Color','r','LineWidth',2)
% % yline(D7_prob_bcell/D7_ss_prod_rate,'Color','r','LineWidth',2)
% plot(true_Bcell_pu1_ss,true_Bcell_pu1_prod,'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
% % text(model_predicted_Bcell_pu1,norm_D7_prob_bcell,"WT PU.1 production in B cell",'FontSize',12,'Color','r','HorizontalAlignment','left')
% 
% % %%%%%% together
% lgd = [lgd,"Remv. rate in B cell","",""];
% % legend(lgd,'location','northwest')
% xlabel('\it{[pu]}')
% ylabel('Normalized Rate')
% set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)
% 
% if contains(plot_option,'log')
%     plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
%     plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
% else
%     plot_ub_x = [0,pu1_totest_range(end)];
%     plot_ub_y = [0,ylim_prod_rate];
% end
% xlim(plot_ub_x)
% ylim(plot_ub_y)
% 
% % title(strrep(fn, '_', ' '))
% % title(final_id)
% % saveas(f,fullfile('outputs_figs',strcat(fn,'_Bcellonly','_',plot_option,'.tif')))
% ax = gca;
% ax.Box = 'on';
% %exportgraphics(f,fullfile('outputs_figs',strcat(fn,'_Bcellonly','_',plot_option,'.tif')),'Resolution',300);

%% plot ProLow 
clear lgd
f = figure('Position',[100 100 800 600]);
axis tight
plot_option = 'log';

% %%%%%%% plot overall patch
patch('Faces',f1,'Vertices',color_patch_pu1_ProLow_overall,'FaceColor',prolowcolor,'EdgeColor',prolowcolor,'FaceAlpha',fa)
hold on
lgd = "Data-constrained pro low";
hold on

% %%%%%%% plot produciont rate for LOW
TR_range_LOW = zeros(size(pu1_totest_range));

for m = 1:length(prolow_all_possibility(:,1))
    cebpa = prolow_all_possibility(m,1);
    % runx1 = prolow_all_possibility(m,2);
    for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
        
        D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        
        D7_curr_prod_rate = D7_prob;
        TR_range_LOW(i) = D7_curr_prod_rate/D7_ss_prod_rate;

    end
    plot(pu1_totest_range,TR_range_LOW,prolow_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
    hold on

    cell_name = prolow_all_possibility_legend{m};
    lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
end

% %%%%%%% generate growth curve
% tau_ProLow = 28.8;
% tau_ProLow_Q1 = 19.9;
% tau_ProLow_Q3 = 47.8;

removal_range_ProLow = rem_rate_constant_ProLow*pu1_totest_range;
min_removal_range_ProLow = rem_rate_constant_ProLow_Q1*pu1_totest_range;
max_removal_range_ProLow = rem_rate_constant_ProLow_Q3*pu1_totest_range;


plot(pu1_totest_range,removal_range_ProLow,'LineWidth',lw,'Color',removalcolor,'LineStyle',removalls)
hold on
plot_range_growthrate_patch(pu1_totest_range,max_removal_range_ProLow,min_removal_range_ProLow,removalcolor)
hold on

% %%%%%% find actual steady state
[~,ind] = mink(abs(TR_range_LOW-removal_range_ProLow),1000);
% concise_ind = extract_intersection_index(ind,30);
true_prolow_pu1_ss = pu1_totest_range(ind(end));
true_prolow_pu1_prod = TR_range_LOW(ind(end));

% %%%%%% mark itnersection
% plot(true_prolow_pu1_ss,true_prolow_pu1_prod,'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
plot(model_predicted_ProLow_pu1,norm_D7_prob_ProLow,'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)

% %%%%%%% together
lgd = [lgd,"Remv. rate in pro low","",""];
legend(lgd,'location','northwest')

xlabel('\it{[pu]}')
ylabel('Normalized Rate')
set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)

if contains(plot_option,'log')
    plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
    plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
else
    plot_ub_x = [0,pu1_totest_range(end)];
    plot_ub_y = [0,ylim_prod_rate];
end
xlim(plot_ub_x)
ylim(plot_ub_y)

% title(strrep(fn, '_', ' '))
% title(final_id)
% saveas(f,fullfile('outputs_figs',strcat(fn,'_ProLow','_',plot_option,'.tif')))
ax = gca;
ax.Box = 'on';
%exportgraphics(f,fullfile('outputs_figs',strcat(fn,'_ProLow','_',plot_option,'.tif')),'Resolution',300);

%% then ProHigh
clear lgd
f = figure('Position',[100 100 800 600]);
axis tight
plot_option = 'linear';

% %%%%%%% plot overall patch
patch('Faces',f1,'Vertices',color_patch_pu1_ProHigh_overall,'FaceColor',prohighcolor,'EdgeColor',prohighcolor,'FaceAlpha',fa)
hold on
lgd = ["Data-constrained pro high"];
hold on

%%%%%%% plot produciont rate for HIGH
TR_range_HIGH = zeros(size(pu1_totest_range));
for m = 1:length(prohigh_all_possibility(:,1))
    cebpa = prohigh_all_possibility(m,1);
    % runx1 = prohigh_all_possibility(m,2);
    for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
        
        D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        
        D7_curr_prod_rate = D7_prob;
        TR_range_HIGH(i) = D7_curr_prod_rate/D7_ss_prod_rate;

    end
    plot(pu1_totest_range,TR_range_HIGH,prohigh_all_possibility_linestyle{m},'LineWidth',lw,'Color',prodcolor)
    hold on

    cell_name = prohigh_all_possibility_legend{m};
    lgd = [lgd,strcat('Prod. rate in'," ",cell_name)];
end


% %%%%%%% generate growth curve
% tau_ProHigh = 11.6;
% tau_ProHigh_Q1 = 9.47;
% tau_ProHigh_Q3 = 15.8;

removal_range_ProHigh = rem_rate_constant_ProHigh*pu1_totest_range;
min_removal_range_ProHigh = rem_rate_constant_ProHigh_Q1*pu1_totest_range;
max_removal_range_ProHigh = rem_rate_constant_ProHigh_Q3*pu1_totest_range;

plot(pu1_totest_range,removal_range_ProHigh,'LineWidth',lw,'Color',removalcolor,'LineStyle',removalls)
hold on
plot_range_growthrate_patch(pu1_totest_range,max_removal_range_ProHigh,min_removal_range_ProHigh,removalcolor)

% %%%%%% find actual steady state
[~,ind] = min(abs(TR_range_HIGH-removal_range_ProHigh));
concise_ind = extract_intersection_index(ind,40);
true_prohigh_pu1_ss = pu1_totest_range(concise_ind);
true_prohigh_pu1_prod = TR_range_HIGH(concise_ind);

% %%%%%% mark itnersection
% plot(true_prohigh_pu1_ss,true_prohigh_pu1_prod,'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
plot(model_predicted_ProHigh_pu1,norm_D7_prob_ProHigh,'LineStyle','none','marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)

% %%%%%%% together
lgd = [lgd,"Remv. rate in pro high","",""];
legend(lgd,'location','southeast')

xlabel('\it{[pu]}')
ylabel('Normalized Rate')
set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)

if contains(plot_option,'log')
    plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
    plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
else
    plot_ub_x = [0,pu1_totest_range(end)];
    plot_ub_y = [0,ylim_prod_rate];
end
xlim(plot_ub_x)
ylim(plot_ub_y)

% title(strrep(fn, '_', ' '))
% title(final_id)
% saveas(f,fullfile('outputs_figs',strcat(fn,'_ProHigh','_',plot_option,'.tif')))
ax = gca;
ax.Box = 'on';
%exportgraphics(f,fullfile('outputs_figs',strcat(fn,'_ProHigh','_',plot_option,'.tif')),'Resolution',300);

%% mac
clear lgd
f = figure('Position',[100 100 800 600]);
axis tight
plot_option = 'linear';

% %%%%%%% plot overall patch
plot(1,1,'LineStyle','none','Marker','s','MarkerSize',ms,'MarkerEdgeColor',maccolor,'MarkerFaceColor',maccolor)
hold on
% lgd = "Data-constrained macrophage";
% hold on

% %%%%%%% plot produciont rate
TR_range = zeros(size(pu1_totest_range));

for m = 1:1
    cebpa = 1;
    % runx1 = 1;
    for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]
        
        D7_prob = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        
        D7_curr_prod_rate = D7_prob;
        TR_range(i) = D7_curr_prod_rate/D7_ss_prod_rate;

    end

    plot(pu1_totest_range,TR_range,'LineWidth',lw,'Color',prodcolor)
    hold on

    cell_name = 'macrophage';
    lgd{m} = strcat('Prod. rate in'," ",cell_name,","," ",'\it[ce]\rm = 1.000');
end


% %%%%%%% generate growth curve
removal_range = rem_rate_constant_Mac*pu1_totest_range;

plot(pu1_totest_range,removal_range,'LineWidth',lw,'Color',removalcolor,'LineStyle',removalls)
hold on
% plot_range_growthrate_patch(pu1_totest_range,max_removal_range,min_removal_range,removalcolor)

% %%%%%%% mark itnersection
plot(1,1,'LineStyle','none','Marker','o','MarkerSize', ms, 'MarkerFaceColor', color_mk_intersect,'MarkerEdgeColor', color_mk_intersect)
% text(1,1,"WT PU.1 production in macrophage",'FontSize',12,'Color','r','HorizontalAlignment','right')

% %%%%%% together
lgd2 = ["",lgd,"Remv. rate in macrophage","",""];
legend(lgd2,'location','northwest')

xlabel('\it{[pu]}')
ylabel('Normalized Rate')
set(gca,'FontSize',fs,'Yscale',plot_option,'Xscale',plot_option)


if contains(plot_option,'log')
    plot_ub_x = [1e-4,10^(log10(pu1_totest_range(end)))];
    plot_ub_y = [1e-4,10^(log10(removal_range(end)+0.5))];
else
    plot_ub_x = [0,pu1_totest_range(end)];
    plot_ub_y = [0,ylim_prod_rate];
end
xlim(plot_ub_x)
ylim(plot_ub_y)

% title(strrep(fn, '_', ' '))
% title(final_id)
% saveas(f,fullfile('outputs_figs',strcat(fn,'_Maconly','_',plot_option,'.tif')))
ax = gca;
ax.Box = 'on';
%exportgraphics(f,fullfile('outputs_figs',strcat(fn,'_Maconly','_',plot_option,'.tif')),'Resolution',300);

%% plot Hill([pu1])/Hill(1) to debug
% figure();
% norm_cebpa_conc = zeros(size(pu1_totest_range));
% for h = 1:length(pu1_totest_range)
%     norm_cebpa_conc(h) = cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_totest_range(h))/cebpa_RegBy_pu1_Hill(p_cebp_hill,1);
% end
% plot(pu1_totest_range,norm_cebpa_conc,'LineWidth',lw)
% xlabel('[pu1]')
% ylabel('Norm cebpa')

%% save best pars
% save(sprintf('cebpa_PU1CEfeedback_prediction_PU1SS_N%.0f.mat',p_cebp_hill(3)), ...
%     'true_Bcell_pu1_ss','true_prolow_pu1_prod','true_prohigh_pu1_prod',...
%     'true_Bcell_pu1_prod','true_prolow_pu1_prod','true_prohigh_pu1_prod','pu1_prod_slope_celltype',...
%     'model_predicted_Bcell_cebpa','model_predicted_ProLow_cebpa','model_predicted_ProHigh_cebpa')
% 
% show_results_real_vs_pred_v2([true_Bcell_pu1_ss,true_prolow_pu1_prod,true_prohigh_pu1_prod],...
%     [true_Bcell_pu1_prod,true_prolow_pu1_prod,true_prohigh_pu1_prod],sprintf('cebpa_PU1CEfeedback_prediction_PU1SS_N%.0f_',p_cebp_hill(3)))