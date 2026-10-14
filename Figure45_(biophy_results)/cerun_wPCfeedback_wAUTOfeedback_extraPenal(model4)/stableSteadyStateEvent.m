function [value,isterminal,direction] = stableSteadyStateEvent(t,y,RemvRatePU1,RemvRateCEBP,p,eq_RUNX)
    
    eq_PU1 = y(1);
    eq_CEBPA = y(2);
    dpudt = fun_dF1dt(RemvRatePU1,p,eq_CEBPA,eq_PU1,eq_RUNX);
    dcedt = fun_dF2dt(RemvRateCEBP,p,eq_PU1,eq_CEBPA);

    fnorm = norm([dpudt,dcedt],Inf);

    stability_score = compute_Jacobian_eig(RemvRatePU1,RemvRateCEBP,p,eq_CEBPA,eq_PU1,eq_RUNX);

    lambda_max = max(real(stability_score));

    tol_f = 1e-8;
    tol_lambda = 1e-2; 

    value = max( ...
        fnorm - tol_f, ...
        lambda_max + tol_lambda );

    isterminal = 1;
    direction = -1; % detecting when real eigenvalue is decreasing from zero

end
