function [x,s,parrange,f] = oneparbif_weigenvalue(odefile,yi,pars,ap,MyContSetOptions)

%set options for equilibrium continuation
opt=contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxStepsize',MyContSetOptions.maxstepsize);
opt=contset(opt,'MaxNumPoints',MyContSetOptions.maxpt);
opt=contset(opt,'MinStepsize',MyContSetOptions.mintepsize);
opt=contset(opt,'Backward',MyContSetOptions.b);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'FunTolerance',10^(-10));
 
%equilibrium continuation
[x0, v0] = init_EP_EP(odefile, yi, pars, ap); % tell MATCONT where to start continuation
[x, v, s, h, f] = cont(@equilibrium, x0, [], opt);
% xb = []; sb = []; % don't need this because x,s,f are arrays already, no need to concatenate 

if isempty(x) || size(x,2) < 2
    warning('Continuation did not converge at the initial point. Skipping this case.');
    parrange = 0;
    % handle gracefully: skip, log, try different x0, etc.
else
    parrange = x(end,:); % all sampled parameter range
    % proceed with normal processing of x, v, s, h, f
end

