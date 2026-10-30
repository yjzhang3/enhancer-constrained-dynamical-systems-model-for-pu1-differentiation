function [final_er,fg_Hill_func] = obj_twoHill_wLG(p)

load(fullfile(fileparts(pwd), 'Error propagation and Normalization', 'error_prop.mat'));

basal = p(1);
vmax = p(2) ;
N = p(3);
khalf = p(4);
g_basal= p(5);
g_vmax=p(6);
g_N = p(7)  ;
g_khalf = p(8);
coop = p(9);
g_input_Bcell = p(10);
g_input_ProLow = p(11);
g_input_ProHigh = p(12);
g_input_mac = p(13);

% function setup
Hill_func = @(x) basal+(vmax-basal)*(x^N/(x^N+khalf^N)); % f(x) is pu1 alone
g_Hill_func = @(w) g_basal+(g_vmax-g_basal)*(w^g_N/(w^g_N+g_khalf^g_N)); % but g(w) should not take 
% the same concentration input as f(x)!! You also can't just use four
% freely optimized constants to multiply onto f(x), because we're assuming
% these four mysterous TF concetration also follows a Hill relationship,
% which is consistent with the biophysical model.
fg_Hill_func = @(x,w) Hill_func(x)*g_Hill_func(w)*coop;

% predict production rate in diff ct
fg_Hill_func_Bcell = fg_Hill_func(conc_ratio_all(1,2),g_input_Bcell);
fg_Hill_func_ProLow = fg_Hill_func(conc_ratio_all(2,2),g_input_ProLow);
fg_Hill_func_ProHigh = fg_Hill_func(conc_ratio_all(3,2),g_input_ProHigh);
fg_Hill_func_mac = fg_Hill_func(1,g_input_mac);
prod_sim_data = [fg_Hill_func_Bcell,fg_Hill_func_ProLow,fg_Hill_func_ProHigh,fg_Hill_func_mac];

% enforce slope of prod < remv (slope with respect to pu1, or x, only. So
% g(w) is a constant to be kept: slope = df(x)/dx*g(w)
Hill_pars = p(1:4);
cebpa_input = p(10:13);
LG_ct = zeros(1,4);
for j = 1:length(LG_ct)
    if j < 4
        LG_ct(j) = Hill_log_gain(Hill_pars,conc_ratio_all(j,2));
    else
        LG_ct(j) = Hill_log_gain(Hill_pars,1);
    end
end

% error (difference between prod and remv)
rem_real_data = [remRate_constant_ratio_all(1,2)*conc_ratio_all(1,2),...
    remRate_constant_ratio_all(2,2)*conc_ratio_all(2,2),...
    remRate_constant_ratio_all(3,2)*conc_ratio_all(3,2),...
    1];
err = 10*weighted_msd(prod_sim_data,rem_real_data); 

% penalty
p1 = penalty_function(g_basal,g_vmax,0.1,0.1);

for i = 1:4
    p1 = p1 + penalty_function(LG_ct(i),1,0.1,0.1);
end

for i = 1:3
    p1 = p1 + penalty_function(prod_sim_data(i),prodRate_ratio_all(i,3),0.1,0.1); % prod < prod upper bound
    p1 = p1 + penalty_function(prodRate_ratio_all(i,1),prod_sim_data(i),0.1,0.1); % prod > prod lower bound
end

final_er = err+p1;

