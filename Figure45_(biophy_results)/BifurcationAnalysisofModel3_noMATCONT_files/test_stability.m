function ind_final = test_stability(par)
% assume par is a n*T matrix of all found fixed points, with n state variables and T variations of
% the control parameter. This function helps find index of all fixed points
% that are stable (all real part of the eigenvalues are negative)

ind_all = zeros(1,length(par(1,:)));
for i = 1:length(ind_all) % each column
    s = 1;
    for j = 1:length(par(:,1)) % look at all state variables (rows)
        s = s*(real(par(j,i))<0);
    end
    ind_all(i) = s;
end
ind_final = find(ind_all); % find all nonzero elements