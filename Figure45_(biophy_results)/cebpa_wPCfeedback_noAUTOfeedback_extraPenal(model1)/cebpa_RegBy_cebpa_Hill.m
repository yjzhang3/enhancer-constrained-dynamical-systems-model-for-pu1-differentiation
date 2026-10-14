function cebpa_conc = cebpa_RegBy_cebpa_Hill(p_auto_hill,ce_conc)
% a phenomnelogical Hill equation that takes in ce concentration and
% outputs cebpa gene expression 
% p_auto_hill is an array containing the parameters of the Hill equation
Amin = p_auto_hill(1);
Amax = p_auto_hill(2);  % we can't let Amax to be 1 because we know maximum
% cepba could be greater than the concentraiton in macrophage (which is
% normalized to 1 now)
N = p_auto_hill(3);
K = p_auto_hill(4);

cebpa_conc = Amin + (Amax-Amin)*ce_conc^N/(K^N+ce_conc^N);

end