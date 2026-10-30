function lg = Hill_log_gain(Hill_pars,x)

p = Hill_pars;
basal = p(1);
vmax = p(2) ;
N = p(3);
khalf = p(4);

output = basal+(vmax-basal)*x^N/(x^N+khalf^N);
slope = (vmax-basal)*N*khalf^N*x^(N-1)/(khalf^N+x^N)^2;
lg = output/x*slope;