function C = C_PU1(C_max,C_min,Kc,N,pu1)

C = C_min + (C_max-C_min)*Kc^N/(Kc^N+pu1^N);