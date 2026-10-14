function [miu,sig] = fit_pars_using_quantile(Qlow,Qhigh,zscore)

miu = (log(Qhigh)+log(Qlow))/2;
sig = (log(Qhigh)-log(Qlow))/(2*zscore);