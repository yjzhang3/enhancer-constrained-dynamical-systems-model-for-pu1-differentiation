%% load no feedback parameters
clear
clc

ct = 2; % pro low
load('best_model_ss.mat','best_pars','runx1_ss_conc','pu1_ss_conc','cebpa_ss_conc')
load('error_prop_0216.mat','remRate_constant_ratio_all','remRate_constant_ratio_all_CEBPA','conc_ratio_all')
energy_table = best_pars(1:16);
nuc_destab_location = [0,1,0,1,0,0];
RNAPHandle = @RNAP_occ_prob_cebpa_runx1_general;
bcellcolor = '#FF5800';
prolowcolor = '#BFBFBF';
prohighcolor = 'k';
color_scheme = {bcellcolor,prolowcolor,prohighcolor};
fa = 0.6; % facealpha of constraint patch (not growth rate patch)
fs = 25;
lg_fs = 18;
lw = 3;
ms = 13;

norm_remv_rate_const = remRate_constant_ratio_all(ct,2);
norm_remv_rate_const_CEBPA = remRate_constant_ratio_all_CEBPA(ct,2);
runx1_cellType_conc = runx1_ss_conc(ct);

%% load mapping between pu and ce based on nofeedback
A = readmatrix("fixedG_1Hill.txt");
[~,I] = min(A(:,end));
p_cebp_hill = A(I,1:end-1);
p_cebp_hill_Amin = p_cebp_hill(1);
p_cebp_hill_Amax = p_cebp_hill(2);
p_cebp_hill_N = p_cebp_hill(3);
p_cebp_hill_K = p_cebp_hill(4);

% load best pars after bifur analysis (now in model 3 I need this, because
% I want to show the best results)
load('CT=2_cebphill_{Amin}_combined.mat')
[minerr,I2] = min(all_err_sorted);
p_cebp_hill(1) = rangealpha_sorted(I2);

%% D7 
site_index = 1;
possible_state_label_14kb = [-14;-1;2;3;14;23;214;314;2314];
possible_state_label_12kb = [-12;-1;12;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_10kb = [-10;-1;10;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_PRO = [-100;-1;100;zeros(length(possible_state_label_14kb)-3,1)];
possible_state_label_persite = [possible_state_label_14kb,possible_state_label_12kb,possible_state_label_10kb,possible_state_label_PRO];
valid_sublist = sublist_finder(possible_state_label_persite,site_index);
D7_valid_sublist = valid_sublist;

%% verify with fsolve
y0 = [0.09,0.09];
ODE_func_ProLow = @(t,y) ODE_func_cellType(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,norm_remv_rate_const_CEBPA,runx1_cellType_conc);
fs_func_ProLow = @(y) ODE_func_cellType([],y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,norm_remv_rate_const_CEBPA,runx1_cellType_conc);
[~,ProLow_conc_over_time] = ode23t(ODE_func_ProLow,[0 400],y0,odeset('RelTol',1e-5,'AbsTol',1e-6));
options = optimset('Display','off');
ProLow_sol = fsolve(fs_func_ProLow,ProLow_conc_over_time(end,:),options);

ProLow_sol_stab = zeros(length(ProLow_sol(:,1)),2);
for i = 1:length(ProLow_sol(:,1))
    ProLow_sol_stab(i,:) = compute_Jacobian_eig_numerical(ODE_func_ProLow,[ProLow_sol(i,1),ProLow_sol(i,2)]);
end
idx2 = find(max(ProLow_sol_stab,[],2) < 0);

%% construct nuclline
close all

% pu1 (y)
dydt = @(x,y) RNAPHandle(D7_valid_sublist,energy_table,[1,y,x,runx1_cellType_conc,1],nuc_destab_location)./...
    RNAPHandle(D7_valid_sublist,energy_table,[1,1,1,1,1],nuc_destab_location) - ... % production rate so far
    norm_remv_rate_const.*y; % and removal rate

% cebpa (x)
dxdt = @(x,y) cebpa_RegBy_pu1_Hill(p_cebp_hill,y)./(cebpa_RegBy_pu1_Hill(p_cebp_hill,1)) ...
    - norm_remv_rate_const_CEBPA.*x;
    % ((log(2)/tau_cellType+D_cebpa)/(log(2)/tau_Mac+D_cebpa))*y(2);
 
f = figure('Position',[100 100 800 600]);
% Plot the x-nullcline (dx/dt = 0)
fhx = fimplicit(@(x,y) dydt(x,y), [0 1], 'Color',[0.0588, 0.2627, 0.5725], 'LineWidth', lw);
fhx_pu = fhx.XData;
fhx_ce = fhx.YData;

hold on
% Plot the y-nullcline (dy/dt = 0)
fhy = fimplicit(@(x,y) dxdt(x,y), [0 1], 'Color',[1.0000, 0.2863, 0.2863], 'LineWidth', lw);
fhy_pu = fhy.XData;
fhy_ce = fhy.YData;
xlabel('\it{[ce]}')
ylabel('\it{[pu]}')

%% overlay experimental [pu]
% patch
hold on;
xl = [0 1];  % get the current x-axis limits [xmin xmax] from the open figure
y0 = conc_ratio_all(ct,1);      % baseline y-value where the patch starts (adjust to your needs)
w = conc_ratio_all(ct,3) - conc_ratio_all(ct,1);
p = patch([xl(1) xl(2) xl(2) xl(1) ],[y0 y0 y0+w y0+w], [0.8 0.8 0.8], ...
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
% clc
intersection_coord = fliplr([0.7571,0.9816]);
ODE_func = @(t,y) ODE_func_cellType(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,norm_remv_rate_const_CEBPA,runx1_cellType_conc);
eigvalue_eval = compute_Jacobian_eig_numerical(ODE_func,[0.7571,0.9816]);
if max(eigvalue_eval) < 0
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',ms*0.8)
else
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','none','MarkerSize',ms*0.8,'LineWidth',3)
end

intersection_coord = fliplr(ProLow_sol);
ODE_func = @(t,y) ODE_func_cellType(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,norm_remv_rate_const_CEBPA,runx1_cellType_conc);
eigvalue_eval = compute_Jacobian_eig_numerical(ODE_func,ProLow_sol);
if max(eigvalue_eval) < 0
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',ms*0.8)
else
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','none','MarkerSize',ms*0.8,'LineWidth',3)
end


hold on
intersection_coord = fliplr([0.2031,0.0917]);
ODE_func = @(t,y) ODE_func_cellType(t,y,RNAPHandle,p_cebp_hill,D7_valid_sublist,energy_table,nuc_destab_location,norm_remv_rate_const,norm_remv_rate_const_CEBPA,runx1_cellType_conc);
eigvalue_eval = compute_Jacobian_eig_numerical(ODE_func,intersection_coord);
if max(eigvalue_eval) < 0
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',ms*0.8)
else
    plot(intersection_coord(1),intersection_coord(2),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','none','MarkerSize',ms*0.8,'LineWidth',3)
end
set(gca,'FontSize',fs)
legend('Measured [pu]','d[pu]/dt=0','d[ce]/dt=0','FontSize',lg_fs,'Location','Northwest')
exportgraphics(f,'Model3_prolow_Nuclline_reversexy.tif','Resolution',300)
saveas(f,'Model3_prolow_Nuclline_reversexy.fig')
%% zoom in
ylim([0.15 0.25])
xlim([0.05 0.15])
legend('off')
exportgraphics(f,'Model3_prolow_Nuclline_Zoomin_reversexy.tif','Resolution',300)
