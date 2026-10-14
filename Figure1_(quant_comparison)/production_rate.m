function Ap = production_rate(p,A_min,A_max,M,Kp)

Ap = A_min + (A_max-A_min)*p.^M./(p.^M+Kp^M);

end