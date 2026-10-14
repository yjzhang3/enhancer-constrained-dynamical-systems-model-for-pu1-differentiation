%% load no feedback parameters
clear
clc

ct = 2; % pro low
load('nofeedback_best_results.mat','best_pars','pu1_ss_conc','cebpa_ss_conc')
load('error_prop_0216.mat','remRate_constant_ratio_all','remRate_constant_ratio_all_CEBPA','conc_ratio_all')
energy_table = best_pars(1:12);
nuc_destab_location = [0,1,0,1];
RNAPHandle = @RNAP_occ_prob_cebpa_general;
norm_remv_rate_const = remRate_constant_ratio_all(ct,2);
norm_remv_rate_const_CEBPA = remRate_constant_ratio_all_CEBPA(ct,2);



bcellcolor = '#FF5800';
prolowcolor = '#BFBFBF';
prohighcolor = 'k';
color_scheme = {bcellcolor,prolowcolor,prohighcolor};
fa = 0.6; % facealpha of constraint patch (not growth rate patch)
fs = 25;
lg_fs = 16;
lw = 3;
ms = 13;

%% load mapping between pu and ce based on nofeedback
A = readmatrix("OneHillPars_fitted_from_noPCfeedback10.txt");
[~,I] = min(A(:,end));
p_cebp_hill = A(I,1:end-1);
p_cebp_hill_Amin = p_cebp_hill(1);
p_cebp_hill_Amax = p_cebp_hill(2);
p_cebp_hill_N = p_cebp_hill(3);
p_cebp_hill_K = p_cebp_hill(4);

% % load best pars after bifur analysis
% load('CT=2_cebphill_{Amin}_combined.mat')
% [minerr,I2] = min(all_err_sorted);
% p_cebp_hill(1) = rangealpha_sorted(I2);

%% D7 
site_index = 1;
possible_state_label_14kb = [-14;-1;2;14;214];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D7_valid_sublist = valid_sublist; % need to define early becuase it'll be used in parallel later

%% verify with fsolve
x0 = [0.1073    0.0901];
ODE_func_Prolow = @(t,y) ODE_func_cellType_noCEAUTO(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,norm_remv_rate_const_CEBPA);
fs_func_Prolow = @(y) ODE_func_cellType_noCEAUTO([],y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,norm_remv_rate_const_CEBPA);
[~,ProLow_conc_over_time] = ode23t(ODE_func_Prolow,[0 400],x0,odeset('RelTol',1e-5,'AbsTol',1e-6));
options = optimset('Display','off');
ProLow_sol = fsolve(fs_func_Prolow,ProLow_conc_over_time(end,:),options);

ProLow_sol_stab = zeros(length(ProLow_sol(:,1)),2);
for i = 1:length(ProLow_sol(:,1))
    ProLow_sol_stab(i,:) = compute_Jacobian_eig_numerical(ODE_func_Prolow,[ProLow_sol(i,1),ProLow_sol(i,2)]);
end
idx2 = find(max(ProLow_sol_stab,[],2) < 0);

%% construct nuclline
close all

% pu1 (x)
dxdt = @(x,y) RNAPHandle(D7_valid_sublist,energy_table,[1,x,y,1],nuc_destab_location)./...
    RNAPHandle(D7_valid_sublist,energy_table,[1,1,1,1],nuc_destab_location) - ... % production rate so far
    norm_remv_rate_const.*x; % and removal rate

% cebpa (y)
dydt = @(x,y) cebpa_RegBy_pu1_Hill(p_cebp_hill,x)./(cebpa_RegBy_pu1_Hill(p_cebp_hill,1)) ...
    - norm_remv_rate_const_CEBPA.*y;
    % ((log(2)/tau_cellType+D_cebpa)/(log(2)/tau_Mac+D_cebpa))*y(2);
 
f = figure('Position',[100 100 800 600]);
% Plot the x-nullcline (dx/dt = 0)
fimplicit(@(x,y) dxdt(x,y), [0 1], 'b', 'LineWidth', 2);
hold on
% Plot the y-nullcline (dy/dt = 0)
fimplicit(@(x,y) dydt(x,y), [0 1], 'r', 'LineWidth', 2);
xlabel('[pu]')
ylabel('[ce]')

%% overlay experimental [pu]
% patch
hold on;
yl = [0 1];  % get the current x-axis limits [xmin xmax] from the open figure
x0 = conc_ratio_all(ct,1);      % baseline y-value where the patch starts (adjust to your needs)
w = conc_ratio_all(ct,3) - conc_ratio_all(ct,1);
p = patch([x0 x0+w x0+w x0], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8], ...
    'FaceAlpha', fa, 'FaceColor',color_scheme{ct},'EdgeColor','none');
uistack(p, 'bottom');   % optional: push patch behind existing data points/lines

%% Add a quiver plot for the phase portrait/vector field
% hold on
% [X, Y] = meshgrid(0:0.1:1, 0:0.1:1); % X and Y matrix should have the same dimension. 
% % Each element in X and each element in Y, AT THE SAME POSITION, combined to become one
% % coordiante. For example, if a=X(1,1) and b=Y(1,1), this corresponds to
% % (a,b) on the Cartesian coordainte with X on horizontal and Y on vertical
% 
% % dx/dt evaluated at (a,b) is the horizontal component of the vector;
% % dy/dt evaluated at (a,b) is the vertical component of the vetor
% U = zeros(size(X));
% V = zeros(size(Y));
% 
% for ii = 1:length(X(:,1))
%     for jj = 1:length(X(1,:))
%         dxdt_eval = dxdt(X(ii,jj),Y(ii,jj));
%         dydt_eval = dydt(X(ii,jj),Y(ii,jj));
%         U(ii,jj) = dxdt_eval;
%         V(ii,jj) = dydt_eval;
%     end
% end
% quiver(X, Y, U, V, 'AutoScale',1,'Color','k');

%% find intersection; compute eigenvalues
% prompt = "What is the intersection? ";
intersection_coord = ProLow_sol;
ODE_func = @(t,y) ODE_func_cellType_noCEAUTO(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,norm_remv_rate_const_CEBPA);
eigvalue_eval = compute_Jacobian_eig_numerical(ODE_func,intersection_coord);
if max(eigvalue_eval) < 0
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',ms*0.8)
else
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','none','MarkerSize',ms*0.8,'LineWidth',3)
end
hold on

intersection_coord = [0.7476,0.9818];
ODE_func = @(t,y) ODE_func_cellType_noCEAUTO(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,norm_remv_rate_const_CEBPA);
eigvalue_eval = compute_Jacobian_eig_numerical(ODE_func,intersection_coord);
if max(eigvalue_eval) < 0
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',ms*0.8)
else
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','none','MarkerSize',ms*0.8,'LineWidth',3)
end
hold on

intersection_coord = [0.203182,0.100934];
ODE_func = @(t,y) ODE_func_cellType_noCEAUTO(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,norm_remv_rate_const_CEBPA);
eigvalue_eval = compute_Jacobian_eig_numerical(ODE_func,intersection_coord);
if max(eigvalue_eval) < 0
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',ms*0.8)
else
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','none','MarkerSize',ms*0.8,'LineWidth',3)
end
set(gca,'FontSize',fs)
legend('Experimental [pu]','d[pu]/dt=0','d[ce]/dt=0','FontSize',lg_fs,'Location','Southeast')
exportgraphics(f,'Model1_prolow_Nuclline.tif','Resolution',300)
saveas(f,'Model1_prolow_Nuclline.fig')

%% zoom in
xlim([0.10 0.25])
ylim([0.05 0.15])
legend('off')
exportgraphics(f,'Model1_prolow_Nuclline_Zoomin.tif','Resolution',300)
