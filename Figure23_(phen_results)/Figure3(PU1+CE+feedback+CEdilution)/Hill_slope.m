function slope = Hill_slope(Hill_pars,x)

p = Hill_pars;
basal = p(1);
vmax = p(2) ;
N = p(3);
khalf = p(4);

slope = (vmax-basal)*N*khalf^N*x^(N-1)/(khalf^N+x^N)^2;