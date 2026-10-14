%%%%%%%%%%%%% in line function
function [err,pred_ce_conc] = find_best_pars(pu1_conc,ce_conc,Hill_pars)
pred_ce_conc = zeros(1,4);
load('error_prop_0216.mat','remRate_constant_ratio_all_CEBPA')
for i = 1:4
    if i == 3
        remRate = remRate_constant_ratio_all_CEBPA(i,1);
    elseif i == 4
        remRate = 1;
    else
        remRate = remRate_constant_ratio_all_CEBPA(i,2);
    end
    pred_ce_conc(i) = cebpa_RegBy_pu1_Hill(Hill_pars,pu1_conc(i))/cebpa_RegBy_pu1_Hill(Hill_pars,pu1_conc(end))/remRate;
end

err = weighted_msd(pred_ce_conc,ce_conc);
end