function p = penalty_function(ksmall,klarge,g1,g2)
% input
% ksmall: parameter that we want to be smaller
% klarge: parameter that we want to be larger
% g1: penalty weight one
% g2: penalty weight two

p = g1*heaviside(ksmall-klarge)*((ksmall-klarge)/klarge)^2 + g2*heaviside(ksmall-klarge);

end