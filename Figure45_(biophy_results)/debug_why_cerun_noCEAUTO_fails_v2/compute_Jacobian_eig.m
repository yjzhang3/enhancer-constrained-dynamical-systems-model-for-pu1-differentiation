function stability_score = compute_Jacobian_eig(RemvRatePU1,RemvRateCEBP,p,eq_CEBPA,eq_PU1,eq_RUNX)

J = [fun_delF1delPU_final(RemvRatePU1,p,eq_CEBPA,eq_PU1,eq_RUNX),fun_delF1delCEBP_final(p,eq_CEBPA,eq_PU1,eq_RUNX);...
    fun_delF2delPU_final(p,eq_PU1,eq_CEBPA),fun_delF2delCEBP_final(RemvRateCEBP,p,eq_PU1,eq_CEBPA)];

stability_score = eig(J);
