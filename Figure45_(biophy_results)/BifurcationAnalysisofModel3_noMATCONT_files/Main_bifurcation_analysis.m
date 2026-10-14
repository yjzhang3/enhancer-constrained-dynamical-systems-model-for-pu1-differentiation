clear
clc
par_labels = { ...
        'cebphill_{Amin}', 'cebphill_{Amax}', 'cebphill_{N}', 'cebphill_{K}', ...
    };

ap_idx_arr = [1,4];
% ap_initial_par_val_arr = [0.1,1,0.1,1];
for ct = 1:1
% for ct = 1:3
    % for ww = 1:length(ap_idx_arr)
    for ww = 2:2
        apn = par_labels{ap_idx_arr(ww)};
        model3_bifur_on_Hill_pars(ct,ap_idx_arr(ww),apn)
    end
end