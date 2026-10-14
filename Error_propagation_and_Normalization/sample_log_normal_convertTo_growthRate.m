function [ratio_mean,ratio_CI] = sample_log_normal_convertTo_growthRate(miu_nonmac,sig_nonmac,miu_mac,sig_mac)

N = 1e6;
nonmac_samples = lognrnd(miu_nonmac,sig_nonmac,[N,1]);
mac_samples = lognrnd(miu_mac,sig_mac,[N,1]);
ratio_samples = (log(2)./nonmac_samples)./(log(2)./mac_samples);

ratio_mean = mean(ratio_samples);
ratio_CI = quantile(ratio_samples, [0.025, 0.975]);