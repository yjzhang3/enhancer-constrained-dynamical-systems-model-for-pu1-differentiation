% error propagation for box plots (PU1 GFP and synthesis rate)
clear
clc
close all
zscore = 0.674; % used to compute standard deviation given the quantile of a normal distribution 
% (XX% quantile is the value below which XX% of the observations lie). In
% other words, the p-th quantile is the inverse of the cumulative probability
% distribution evaluated at p (i.e., the value at which the cumulative
% probability is less than p). Each quantile has a corresponding Z score,
% which tells you how far away a value is from the mean of the standard
% deviation.

%% data reported range; each arr = [25% quantile, median, 75% quantile]

% PU1-GFP
% pu1_bcell_conc_range = [1/200,0.84/7.42]; 
bcell_conc_range = [0.50,0.62,0.84]*10^4;
ProLow_conc_range = [1.35,1.63,2.36]*10^4; 
ProHigh_conc_range = [3.53,3.78,3.99]*10^4;
nonmac_conc_range = [bcell_conc_range;ProLow_conc_range;ProHigh_conc_range];
mac_conc_range = [6.92,7.42,8.46]*10^4;

% synthesis/production rate of PU1
bcell_prod_rate_range = [0.051,0.10,0.18]*10^3;
ProLow_prod_rate_range = [0.31,0.44,0.69]*10^3;
ProHigh_prod_rate_range = [2.06,2.54,2.72]*10^3;
nonmac_prod_range = [bcell_prod_rate_range;ProLow_prod_rate_range;ProHigh_prod_rate_range];
mac_prod_rate_range = [1.15,2.30,3.33]*10^3;

% assume the *log* of the P(X) (PDF for PU1-GFP or synthesis rate) is nomrally
% distributed, which is consistent with the skewed distribution
% miu of log(P(X)) = median of log(P(X)) = log of median of P(X)
% sigma of log(P(X)) = ( log(Q75)-log(Q25) ) / (2*0.674), where QXX is the
% quantile of P(X)

%% error prop for PU1GFP
conc_ratio_all = zeros(3,3);
cellName_all = {'B cell','pro low','pro high'};
ylabelName = 'PU.1-GFP (a.u.)';
for ct = 1:3
    nonmac_arr = nonmac_conc_range(ct,:);
    mac_arr = mac_conc_range;
    [miu_nonmac,sig_nonmac] = fit_pars_using_quantile(nonmac_arr(1),nonmac_arr(3),zscore);
    [miu_mac,sig_mac] = fit_pars_using_quantile(mac_arr(1),mac_arr(3),zscore);
    [ratio_median,ratio_2575] = sample_log_normal(miu_nonmac,sig_nonmac,miu_mac,sig_mac);
    sample_log_normal_NoRatio_boxPlot(miu_nonmac,sig_nonmac,miu_mac,sig_mac,cellName_all{ct},ylabelName)
    conc_ratio_all(ct,2) = ratio_median;
    conc_ratio_all(ct,[1,3]) = ratio_2575;
end

%% error prop for synthesis rate
close all
prodRate_ratio_all = zeros(3,3);
cellName_all = {'B cell','pro low','pro high'};
ylabelName = 'PU.1 synthesis rate (a.u. hr^{-1})';
for ct = 1:3
    nonmac_arr = nonmac_prod_range(ct,:);
    mac_arr = mac_prod_rate_range;
    [miu_nonmac,sig_nonmac] = fit_pars_using_quantile(nonmac_arr(1),nonmac_arr(3),zscore);
    [miu_mac,sig_mac] = fit_pars_using_quantile(mac_arr(1),mac_arr(3),zscore);
    [ratio_median,ratio_2575] = sample_log_normal(miu_nonmac,sig_nonmac,miu_mac,sig_mac);
    sample_log_normal_NoRatio_boxPlot(miu_nonmac,sig_nonmac,miu_mac,sig_mac,cellName_all{ct},ylabelName)
    prodRate_ratio_all(ct,2) = ratio_median;
    prodRate_ratio_all(ct,[1,3]) = ratio_2575;
end

%% save
conc_ratio_all(1,1) = 1/200;
save('error_prop.mat','conc_ratio_all','prodRate_ratio_all');