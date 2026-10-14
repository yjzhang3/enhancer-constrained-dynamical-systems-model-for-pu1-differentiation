function cebpa_prodRate = cebpa_RegBy_pu1_Hill(p_cebp_hill,pu1_conc)
% a phenomnelogical Hill equation that takes in pu1 concentration and
% outputs cebpa gene expression (assuming to be proportional to C/EBPa
% protein concentration that will affect PU1 gene expression, completing
% the positive feedback loop)
% p_cebp_hill is an array containing the parameters of the Hill equation
Amin = p_cebp_hill(1);
Amax = p_cebp_hill(2);  % we can't let Amax to be 1 because we know maximum
% cepba could be greater than the concentraiton in macrophage (which is
% normalized to 1 now)
N = p_cebp_hill(3);
K = p_cebp_hill(4);

cebpa_prodRate = Amin + (Amax-Amin)*pu1_conc^N/(K^N+pu1_conc^N);

end