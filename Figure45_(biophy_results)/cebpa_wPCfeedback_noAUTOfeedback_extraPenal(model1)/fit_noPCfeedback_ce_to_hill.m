function [err,pred_ce_conc] = fit_noPCfeedback_ce_to_hill(p,pu1_conc,real_ce_conc)
% real_ce_conc are the ones fitted from no feedback model. THis function is
% used to examine the phen. property of hill equation such that stable
% steady state of pu1 can be predicted accurately.

Amin = p(1);
Amax = p(2);
N = p(3);
K = p(4);

Hill_output = @(Amin,Amax,N,K,pu1_conc) Amin + (Amax-Amin)*pu1_conc^N/(K^N+pu1_conc^N);

pu1_conc_bcell = pu1_conc(1);
pu1_conc_ProLow = pu1_conc(2);
pu1_conc_ProHigh = pu1_conc(3);

pred_ce_conc_bcell = Hill_output(Amin,Amax,N,K,pu1_conc_bcell);
pred_ce_conc_ProLow = Hill_output(Amin,Amax,N,K,pu1_conc_ProLow);
pred_ce_conc_ProHigh = Hill_output(Amin,Amax,N,K,pu1_conc_ProHigh);

pred_ce_conc = [pred_ce_conc_bcell,pred_ce_conc_ProLow,pred_ce_conc_ProHigh];

err = weighted_msd(pred_ce_conc/Hill_output(Amin,Amax,N,K,1),real_ce_conc);