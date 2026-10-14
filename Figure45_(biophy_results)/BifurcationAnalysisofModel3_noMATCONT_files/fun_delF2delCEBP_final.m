function final = fun_delF2delCEBP_final(RemvRateCEBP,p,eq_PU1,eq_CEBPA)

load('delF2.mat','fun_F2normconstant','func_delF2delCEBP_raw')
p_cebp_hill = p(20:23);
p_auto_hill = p(24:27);

Amin1 = p_cebp_hill(1);
Amax1 = p_cebp_hill(2);
N1 = p_cebp_hill(3);
K1 = p_cebp_hill(4);

Amin2 = p_auto_hill(1);
Amax2 = p_auto_hill(2);
N2 = p_auto_hill(3);
K2 = p_auto_hill(4);

raw = func_delF2delCEBP_raw(Amax1,Amax2,Amin1,Amin2,eq_CEBPA,K1,K2,N1,N2,eq_PU1);
norm_constant = fun_F2normconstant(Amax1,Amax2,Amin1,Amin2,1,K1,K2,N1,N2,1);
final = (1/norm_constant)*raw-RemvRateCEBP;