function [final_er,fg_Hill_func] = obj_twoHill(p)

load(fullfile(fileparts(pwd), 'Figure1A', 'error_prop.mat'));

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
g_Hill_func = @(w) g_basal+(g_vmax-g_basal)*(w^g_N/(w^g_N+g_khalf^g_N)); % but g(x) should not take 
% the same concentration input as f(x)!!
fg_Hill_func = @(x,w) Hill_func(x)*g_Hill_func(w)*coop;

% predict production rate in diff ct
fg_Hill_func_Bcell = fg_Hill_func(conc_ratio_all(1,2),g_input_Bcell);
fg_Hill_func_ProLow = fg_Hill_func(conc_ratio_all(2,2),g_input_ProLow);
fg_Hill_func_ProHigh = fg_Hill_func(conc_ratio_all(3,2),g_input_ProHigh);
fg_Hill_func_mac = fg_Hill_func(1,g_input_mac);
prod_sim_data = [fg_Hill_func_Bcell,fg_Hill_func_ProLow,fg_Hill_func_ProHigh,fg_Hill_func_mac];

% error
grow_real_data = [remRate_constant_ratio_all(1,2)*conc_ratio_all(1,2),...
    remRate_constant_ratio_all(2,2)*conc_ratio_all(2,2),...
    remRate_constant_ratio_all(3,2)*conc_ratio_all(3,2),...
    1];
prod_real_data = [prodRate_ratio_all(1,2),prodRate_ratio_all(2,2),prodRate_ratio_all(3,2),1];

err = weighted_msd(prod_sim_data,grow_real_data)+weighted_msd(prod_sim_data,prod_real_data);

% penalty
p1 = penalty_function(g_basal,g_vmax,0.1,0.1);

final_er = err+p1;
end

