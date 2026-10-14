function [final_err] = obj_ps_find_Hill_pars(p)

hill_pars = p(1:4);
ct_pu = p(5:8);
load('error_prop_0216.mat');

% hill_pars = basal,vamx, N, k

% production rate
Hill_func = @(hill_pars,x) hill_pars(1)+(hill_pars(2)-hill_pars(1))*(x^hill_pars(3)/(x^hill_pars(3)+hill_pars(4)^hill_pars(3)));
pred_prod_Bcell = Hill_func(hill_pars,ct_pu(1));
pred_prod_ProLow = Hill_func(hill_pars,ct_pu(2));
pred_prod_ProHigh = Hill_func(hill_pars,ct_pu(3));
pred_prod_mac = Hill_func(hill_pars,ct_pu(4));

% slope of a Hill function in terms of x: (v-b)*k^n*n*x^(n-1)/(k^n+x^n)^2
Hill_slope_func = @(hill_pars,x) (hill_pars(2)-hill_pars(1))*hill_pars(4)^hill_pars(3)*hill_pars(3)*x^(hill_pars(3)-1)/...
    (hill_pars(4)^hill_pars(3)+x^hill_pars(3))^2;
slope_pred_prod_Bcell = Hill_slope_func(hill_pars,ct_pu(1));
slope_pred_prod_ProLow = Hill_slope_func(hill_pars,ct_pu(2));
slope_pred_prod_ProHigh = Hill_slope_func(hill_pars,ct_pu(3));
slope_pred_prod_mac = Hill_slope_func(hill_pars,ct_pu(4));

% removal rate
remv_Bcell = remRate_constant_ratio_all(1,2)*ct_pu(1);
remv_ProLow = remRate_constant_ratio_all(2,2)*ct_pu(2);
remv_ProHigh = remRate_constant_ratio_all(3,1)*ct_pu(3);
remv_Mac = 1*ct_pu(4);

% prod = remov
er1 = weighted_msd([pred_prod_Bcell,pred_prod_ProLow,pred_prod_ProHigh,pred_prod_mac],...
    [remv_Bcell,remv_ProLow,remv_ProHigh,remv_Mac]);

% penalty: penalty_function(ksmall,klarge,g1,g2)
g1 = 0.1; g2 = 0.1;
p1 = penalty_function(prodRate_ratio_all(1,1),pred_prod_Bcell,g1,g2) + ...
    penalty_function(pred_prod_Bcell,prodRate_ratio_all(1,3),g1,g2);
p2 = penalty_function(prodRate_ratio_all(2,1),pred_prod_ProLow,g1,g2) + ...
    penalty_function(pred_prod_ProLow,prodRate_ratio_all(2,3),g1,g2);
p3 = penalty_function(prodRate_ratio_all(3,1),pred_prod_ProHigh,g1,g2) + ...
    penalty_function(pred_prod_ProHigh,prodRate_ratio_all(3,3),g1,g2);
p4 = penalty_function(slope_pred_prod_Bcell,remRate_constant_ratio_all(1,2),g1,g2);
p5 = penalty_function(slope_pred_prod_ProLow,remRate_constant_ratio_all(2,2),g1,g2);
p6 = penalty_function(slope_pred_prod_ProHigh,remRate_constant_ratio_all(3,1),g1,g2);

% together
final_err = er1 + p1 + p2 + p3 + p4 + p5 + p6;