%% D7 all possible configurations
clear
clc

site_index = 1;
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D7_valid_sublist = valid_sublist;

% expression for all configurations
expr_list = cell(length(D7_valid_sublist(:,1)),1);
nuc_destab_location = [0,1,0,1,0,0];
for nn = 1:length(D7_valid_sublist(:,1))
    final_expr_str = BF_per_state_with_cebpa_runx1_general_expression(D7_valid_sublist(nn,:),nuc_destab_location);
    expr_list{nn,1} = final_expr_str;
end

% expression RNAP bound configurations 
rows_with_100 = find(any(D7_valid_sublist == 100, 2));  % logical vector
RNAPbound_expr_list = expr_list(rows_with_100);

%% partition function expression
Z = strcat('(',strjoin(expr_list,'+'),')');
numerator = strcat('(',strjoin(RNAPbound_expr_list,'+'),')');
prob_expr = strcat(numerator,'/',Z);

%% convert to symbolic
eqn = str2sym(prob_expr);
var = symvar(eqn);

delF1delPU_raw = diff(eqn,var(18)); % upper left (need to multiply this by the dimensionless constant, and then subtract the removal rate  of PU)
delF1delCEBP_raw = diff(eqn,var(1)); % upper right (need to multiply this by the dimensionless constant)

fun_F1normconstant = matlabFunction(eqn); % this is the probability expression byitself. It 
% can be used to compute either the normalization constant, or the actual
% probability
fun_delF1delPU_raw = matlabFunction(delF1delPU_raw);
fun_delF1delCEBP_raw = matlabFunction(delF1delCEBP_raw);

save('delF1','fun_F1normconstant','fun_delF1delCEBP_raw','fun_delF1delPU_raw')
% see fun_delF1delPU_final for the final results

%% What about delF2 del?
CE_PROD_expr = '(Amin1 + (Amax1 - Amin1)*PU1^N1/(PU1^N1+K1^N1))*(Amin2 + (Amax2 - Amin2)*CEBPA^N2/(CEBPA^N2+K2^N2))';
eqn = str2sym(CE_PROD_expr);
var  = symvar(eqn);
delF2delPU_raw = diff(eqn,var(end)); % lower left (need to multiply by the dimensionless constant)
delF2delCEBP_raw = diff(eqn,var(5));  % lower right (need to multiply by the dimensionless constant, and subtract the removal rate of cebpa)

fun_F2normconstant = matlabFunction(eqn);
func_delF2delPU_raw = matlabFunction(delF2delPU_raw);
func_delF2delCEBP_raw = matlabFunction(delF2delCEBP_raw);
save('delF2','func_delF2delPU_raw','func_delF2delCEBP_raw','fun_F2normconstant')