%%%%%%%%%%%%%%%%%%%%% debug why CEAUTO is necessary

%% using best fitted pars to see if Hill([pu]) is enough
clear
clc
% delete HillPars_fitted_from_model4.txt

% load best pars from model 4
% load('best_model_ss.mat')
load('best_model_ss_rescaLS1.mat')

% find a set of Hill pars that satisfy best pu1_ss map to best ce_ss
objFunc = @(Hill_pars) find_best_pars(pu1_ss_conc,cebpa_ss_conc,Hill_pars);

% use the same bound for model 3
lb = [0.0001,1,2,0.001];
ub = [0.1,1,2,2];
% fn = 'HillPars_fitted_from_model4.txt';
fn = 'HillPars_fitted_from_model4_rescaleLS1.txt';

for i = 1:10
[x,fval] = particleswarm(objFunc,4,lb,ub);
create_parameter_file(fn,x,fval,'a+')
end

% extract best pars
A = readmatrix(fn);
[minerr,I] = min(A(:,end));
best_pars_Hill = A(I,1:end-1);
[err,pred_ce_conc] = find_best_pars(pu1_ss_conc,cebpa_ss_conc,best_pars_Hill);

% plot
f=figure();
plot(pu1_ss_conc,cebpa_ss_conc,'LineStyle','none','Marker','o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',15)
hold on
plot(pu1_ss_conc,pred_ce_conc,'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',15)
legend('Best-fitted [ce]','Predicted [ce], N = 2','Location','best')
xlabel('Best-fitted [pu]')
ylabel('[ce]')
set(gca,'FontSize',20)
exportgraphics(f,'OneHillPars_fitted_from_model4_bad.tif','Resolution',300)

%% ok, so far we see that N = 2 of Hill([pu]) can't capture the required
% trend of [ce]

% what if I find the good pars, and plug that back into the original
% objective function. Will stability change? Otherwise, a higher N of
% Hill([pu]) is also ok?

clear
clc
% delete HillPars_fitted_from_model4.txt

% load best pars from model 4
% load('best_model_ss.mat')
load('best_model_ss_rescaLS1.mat')


% find a set of Hill pars that satisfy best pu1_ss map to best ce_ss
objFunc = @(Hill_pars) find_best_pars(pu1_ss_conc,cebpa_ss_conc,Hill_pars);

% use the same bound for model 3
lb = [0.0001,1,15,0.001];
ub = [0.1,1,15,2];
% fn = 'HillPars_fitted_from_model4_Nfree.txt';
fn = 'HillPars_fitted_from_model4_Nfree_rescaleLS1.txt';

for i = 1:10
[x,fval] = particleswarm(objFunc,4,lb,ub);
create_parameter_file(fn,x,fval,'a+')
end

% extract best pars
A = readmatrix(fn);
[minerr,I] = min(A(:,end));
best_pars_Hill = A(I,1:end-1);
[err,pred_ce_conc] = find_best_pars(pu1_ss_conc,cebpa_ss_conc,best_pars_Hill);

% plot
f=figure();
plot(pu1_ss_conc,cebpa_ss_conc,'LineStyle','none','Marker','o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',15)
hold on
plot(pu1_ss_conc,pred_ce_conc,'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',15)
legend('Best-fitted [ce]',strcat('Predicted [ce]',',',sprintf('N = %.3f',best_pars_Hill(3))),'location','best')
xlabel('Best-fitted [pu]')
ylabel('[ce]')
set(gca,'FontSize',20)
exportgraphics(f,'OneHillPars_fitted_from_model4.tif','Resolution',300)


%% first test whether this can replicate the error on PU.1 ss, PU.1 prod raate, and gene expression data
% load best pars from model 4
load('best_model_ss.mat')
debug_no_ODE_test_error(pu1_ss_conc,pred_ce_conc,runx1_ss_conc,energy_table,best_pars(17))

%% now test stability using the original objective function
% using best-fitted pu1 and cebpa as initial condition, and use independtly
% fitted Hill([pu]) pars, see if the predicteed ce is stable
[tb,Bcell_conc_over_time,tpl,ProLow_conc_over_time,tph,ProHigh_conc_over_time] = ...
    debug_ODE_stability(pu1_ss_conc,cebpa_ss_conc,best_pars_Hill,energy_table,runx1_ss_conc);

% color for constraint patch
bcellcolor = '#FF5800';
prolowcolor = '#BFBFBF';
prohighcolor = 'k';
maccolor = '#2626B2';

f=figure('position',[148,278,1088,420]);
subplot(1,2,1)
plot(tb,Bcell_conc_over_time(:,2),'LineWidth',3,'Color',bcellcolor)
hold on
plot(tpl,ProLow_conc_over_time(:,2),'LineWidth',3,'Color',prolowcolor)
hold on
plot(tph,ProHigh_conc_over_time(:,2),'LineWidth',3,'Color',prohighcolor)
legend('Bcell','ProLow','ProHigh')
xlabel('time')
ylabel('d[ce]/dt')
set(gca,'FontSize',20)

subplot(1,2,2)
plot(tb,Bcell_conc_over_time(:,1),'LineWidth',3,'Color',bcellcolor)
hold on
plot(tpl,ProLow_conc_over_time(:,1),'LineWidth',3,'Color',prolowcolor)
hold on
plot(tph,ProHigh_conc_over_time(:,1),'LineWidth',3,'Color',prohighcolor)
legend('Bcell','ProLow','ProHigh')
xlabel('time')
ylabel('d[pu]/dt')
set(gca,'FontSize',20)
exportgraphics(f,'Plug_OneHillPars_fitted_from_model4.tif','Resolution',300)
