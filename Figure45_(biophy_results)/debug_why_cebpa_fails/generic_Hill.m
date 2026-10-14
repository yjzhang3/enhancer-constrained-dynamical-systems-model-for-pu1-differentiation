function cebpa_conc = generic_Hill(p_hill,tf_conc)
% a phenomnelogical Hill equation that takes in ce concentration and
% outputs cebpa gene expression 
% p_auto_hill is an array containing the parameters of the Hill equation
Amin = p_hill(1);
Amax = p_hill(2);  % we can't let Amax to be 1 because we know maximum
% cepba could be greater than the concentraiton in macrophage (which is
% normalized to 1 now)
N = p_hill(3);
K = p_hill(4);

cebpa_conc = Amin + (Amax-Amin)*tf_conc^N/(K^N+tf_conc^N);

end