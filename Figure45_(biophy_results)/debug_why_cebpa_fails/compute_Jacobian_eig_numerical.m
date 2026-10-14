function stability_score = compute_Jacobian_eig_numerical(FunHandle,y_eq)
J = numericalJacobian(@(y) FunHandle([],y), y_eq);
stability_score = eig(J);
end

function J = numericalJacobian(fun,x)
% chatGPT's recommendation to numerically 

n = length(x);
J = zeros(n);
h = 1e-6;

for i = 1:n
    xp = x;
    xm = x;

    xp(i) = xp(i) + h;
    xm(i) = xm(i) - h;

    J(:,i) = (fun(xp) - fun(xm))/(2*h);
end

end