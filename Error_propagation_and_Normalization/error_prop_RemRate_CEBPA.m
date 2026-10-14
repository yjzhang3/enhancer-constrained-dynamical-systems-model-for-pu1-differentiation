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

remRate_constant_ratio_all_CEBPA = zeros(3,3);

for ct = 1:3
    nonmac_arr = nonmac_tau(ct,:);
    mac_arr = tau_Mac;
    [miu_nonmac,sig_nonmac] = fit_pars_using_quantile(nonmac_arr(1),nonmac_arr(3),zscore);
    [miu_mac,sig_mac] = fit_pars_using_quantile(mac_arr(1),mac_arr(3),zscore);
    [ratio_mean_rate,ratio_CI_rate] = sample_log_normal_convertTo_RemRate_CEBPA(miu_nonmac,sig_nonmac,miu_mac,sig_mac);
    remRate_constant_ratio_all_CEBPA(ct,2) = ratio_mean_rate;
    remRate_constant_ratio_all_CEBPA(ct,[1,3]) = ratio_CI_rate;
end

save('error_prop_0216.mat','remRate_constant_ratio_all_CEBPA','-append')