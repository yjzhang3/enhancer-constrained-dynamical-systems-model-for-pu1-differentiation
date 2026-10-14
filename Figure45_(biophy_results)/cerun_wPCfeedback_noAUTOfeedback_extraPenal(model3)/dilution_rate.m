function Cp = dilution_rate(p,C_max,C_min,N,Kc)

Cp = C_min + (C_max-C_min)*Kc^N/(Kc^N+p^N);

end