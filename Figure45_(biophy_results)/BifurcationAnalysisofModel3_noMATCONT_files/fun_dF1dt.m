function final = fun_dF1dt(RemvRatePU1,p,eq_CEBP,eq_PU1,eq_RUNX) 
% only applicable to objective_2D_cerun_PCfeed_CEauto
load('delF1.mat','fun_F1normconstant')

G_NUC = p(1);
G_PU1_14 = p(2);
G_PU1_12 = p(3);
G_PU1_10 = p(4);
G_CEBP = p(5);
G_RUNX = p(6);
G_RNAP = p(7);
G_NUC_PU1 = p(8); 
G_NUC_CEBP = p(9);
G_NUC_RUNX1 = p(10);
G_PU1_RNAP = p(11);
G_CEBP_RNAP = p(12);
G_RUNX1_RNAP = p(13);
G_RUNX_PU1 = p(14);
G_RUNX_CEBP = p(15);
G_PU1_CEBP = p(16);

norm_constant = fun_F1normconstant(1,G_CEBP,G_CEBP_RNAP,G_NUC,G_NUC_PU1,G_NUC_CEBP,G_PU1_10,G_PU1_12,G_PU1_14,G_PU1_CEBP,G_PU1_RNAP,G_RNAP,G_RUNX,G_RUNX_PU1,G_RUNX_CEBP,G_RUNX1_RNAP,1,1,1,1);
raw = fun_F1normconstant(eq_CEBP,G_CEBP,G_CEBP_RNAP,G_NUC,G_NUC_PU1,G_NUC_CEBP,G_PU1_10,G_PU1_12,G_PU1_14,G_PU1_CEBP,G_PU1_RNAP,G_RNAP,G_RUNX,G_RUNX_PU1,G_RUNX_CEBP,G_RUNX1_RNAP,1,eq_PU1,1,eq_RUNX);

final = (1/norm_constant)*raw - RemvRatePU1*eq_PU1;