function [ratio_median,ratio_CI] = sample_log_normal(miu_nonmac,sig_nonmac,miu_mac,sig_mac)

N = 1e6;
nonmac_samples = lognrnd(miu_nonmac,sig_nonmac,[N,1]);
mac_samples = lognrnd(miu_mac,sig_mac,[N,1]);
ratio_samples = nonmac_samples./mac_samples;

ratio_median = median(ratio_samples);
ratio_CI = quantile(ratio_samples, [0.25, 0.75]);