% error propagation in cell cycle length
clear
clc
zscore = 1.96;

%% reported measurement; arr = [lower 95CI, mean, higher 95CI]
tau_Mac = [37.8,40,56.4];
tau_Bcell = [17.7,25.5,41.2];
tau_ProLow = [19.9,28.8,47.8];
tau_ProHigh = [9.47,11.6,15.8];
nonmac_tau = [tau_Bcell;tau_ProLow;tau_ProHigh];

%% error prop for cell cycle
growthRate_ratio_all = zeros(3,3);
remRate_constant_ratio_all = zeros(3,3);

cellcycle_length_ratio_all = zeros(3,3);
cellName_all = {'B cell','pro low','pro high'};
ylabelName = 'Cell cycle length (hr)';

for ct = 1:3
    nonmac_arr = nonmac_tau(ct,:);
    mac_arr = tau_Mac;
    [miu_nonmac,sig_nonmac] = fit_pars_using_quantile(nonmac_arr(1),nonmac_arr(3),zscore);
    [miu_mac,sig_mac] = fit_pars_using_quantile(mac_arr(1),mac_arr(3),zscore);
    sample_log_normal_NoRatio_95CI(miu_nonmac,sig_nonmac,miu_mac,sig_mac,cellName_all{ct},ylabelName)

    [ratio_mean_length,ratio_CI_length] = sample_log_normal(miu_nonmac,sig_nonmac,miu_mac,sig_mac);
    cellcycle_length_ratio_all(ct,2) = ratio_mean_length;
    cellcycle_length_ratio_all(ct,[1,3]) = ratio_CI_length;

    [ratio_mean_rate,ratio_CI_rate] = sample_log_normal_convertTo_growthRate(miu_nonmac,sig_nonmac,miu_mac,sig_mac);
    % the ratio here is directly growth rate of non-mac/growth rate of mac,
    % not the actual normalized removal rate in the dimensionless ODE (so
    % can't be used to plot in results!!). These are only for fig1
    growthRate_ratio_all(ct,2) = ratio_mean_rate;
    growthRate_ratio_all(ct,[1,3]) = ratio_CI_rate;
    
    [ratio_mean_rate,ratio_CI_rate] = sample_log_normal_convertTo_RemRate(miu_nonmac,sig_nonmac,miu_mac,sig_mac);
    remRate_constant_ratio_all(ct,2) = ratio_mean_rate;
    remRate_constant_ratio_all(ct,[1,3]) = ratio_CI_rate;
end

%% save
save('error_prop.mat','growthRate_ratio_all','remRate_constant_ratio_all','cellcycle_length_ratio_all','-append')

%% compute the ratios between removal rates across cell types
remRate_PL_to_Bcell = remRate_constant_ratio_all(2,2)/remRate_constant_ratio_all(1,2);
remRate_PH_to_PL = remRate_constant_ratio_all(3,2)/remRate_constant_ratio_all(2,2);
remRate_mac_to_PH = 1/remRate_constant_ratio_all(3,2);

% the minimum PU.1 production rate of each cell type for steady-state PU.1 to have an
% increasing level from B cell to mac
min_PH_prodRate = 1/remRate_mac_to_PH; % the production rate in mac should be 1
min_PL_prodRate = min_PH_prodRate/remRate_PH_to_PL;
min_Bcell_prodRate = min_PL_prodRate/remRate_PL_to_Bcell;

% ok, even if min_PL_prodRate < min_Bcell_prodRate
% but also min_PL_prodRate < min_PH_prodRate, indicating that we don't
% need a monotonic increasing production rate, because we said PU.1
% production rate is Hill([pu1])*Hill([X]), and we want the entire
% production rate is monotonic with respect to [pu1], then [X] must also be
% monotonic 