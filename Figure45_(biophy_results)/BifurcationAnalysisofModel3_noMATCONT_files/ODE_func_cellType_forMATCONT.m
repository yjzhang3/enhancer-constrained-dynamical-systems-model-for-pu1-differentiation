function dydt = ODE_func_cellType_forMATCONT(t,y,varargin)

%% pass in extra parameters exactly as defined in init_EP
% order of parameters in your pars vector must exactly match the order you deal them into p1, p2, ... — MATCONT has no idea what your parameters are named, it only passes them positionally.
pars = varargin;
[p_cebp_hill_Amin,p_cebp_hill_Amax,p_cebp_hill_N,p_cebp_hill_K,ct] = deal(pars{:});

p_cebp_hill = [p_cebp_hill_Amin,p_cebp_hill_Amax,p_cebp_hill_N,p_cebp_hill_K];

% p_auto_hill = [p_auto_hill_Amin,p_auto_hill_Amax,p_auto_hill_N,p_auto_hill_K];

%% RNAP handle
RNAPHandle = @RNAP_occ_prob_cebpa_runx1_general;

%% D7 configurations
site_index = 1;
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D7_valid_sublist = valid_sublist;

%% nuc_destab_location
nuc_destab_location = [0,1,0,1,0,0];

%% norm_remv_rate_const
load('error_prop_0216.mat','remRate_constant_ratio_all','remRate_constant_ratio_all_CEBPA')
norm_remv_rate_const = remRate_constant_ratio_all(ct,2);
norm_remv_rate_const_CEBPA = remRate_constant_ratio_all_CEBPA(ct,2);

if ct == 3
    norm_remv_rate_const = remRate_constant_ratio_all(ct,1);
    norm_remv_rate_const_CEBPA = remRate_constant_ratio_all_CEBPA(ct,1);
end
%% runx1 concentration and energy table
load("best_model_ss.mat",'energy_table','runx1_ss_conc')
runx1_cellType_conc = runx1_ss_conc(ct);

% pu1
dpudt = RNAPHandle(D7_valid_sublist,energy_table,[1,y(1),y(2),runx1_cellType_conc,1],nuc_destab_location)/...
    RNAPHandle(D7_valid_sublist,energy_table,[1,1,1,1,1],nuc_destab_location) - ... % production rate so far
    norm_remv_rate_const*y(1); % and removal rate

% cebpa
dcedt = cebpa_RegBy_pu1_Hill(p_cebp_hill,y(1))/(cebpa_RegBy_pu1_Hill(p_cebp_hill,1)) ...
    - norm_remv_rate_const_CEBPA*y(2); % production rate so far
%     ((log(2)/tau_cellType+D_cebpa)/(log(2)/tau_Mac+D_cebpa))*y(2);

dydt = [dpudt;dcedt];

