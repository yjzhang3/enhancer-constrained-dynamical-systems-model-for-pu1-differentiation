% v2 is using yesfeedbacknoODE results for initial condition (that is, from
% local_fit_cebpa_general_PU1CEfeedback_smallestPHdilute_0401_PattrC_freeN_maxV1_1234567891011_0101)

% v3 is the same as v2, except I break simulation into two separate times

% v4 is using true steady state of v3 (that is, manulally add more time to
% v3 best results to find the true steady state) as initial condition; also
% break simulation into two separate times, with the second simulation much
% longer (most imporatntly, increased tolerance!!)

% v5 has same initial conditions as v3, and I increase tolerance for both
% simulations