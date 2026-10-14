function ODE_fit_with_cebpa_runx1_PU1CEfeedback_GA(objFuncHandle,nuc_destab_location,selec)

%% housekeeping
load(fullfile(parentdir(pwd, 1), 'Error propagation and Normalization','error_prop.mat')); 
fun = @(p) objFuncHandle(selec,p,nuc_destab_location,conc_ratio_all,prodRate_ratio_all);
funname = func2str(objFuncHandle);

orig_name = strcat('ODElocalfit_GA_',funname(15:end),'_PrepC','_2N','_maxV1'); 

fn = strcat(orig_name,'_',num2str(selec),'_',num2str(nuc_destab_location));
fn = strrep(fn, ' ', '');
final_fn = strcat(fn,'.txt');
% task_fn = fullfile(NOTS_path,final_fn);
nvars = 24;

%% parameter bound
% use previous model results to help optimize faster
% previous_results = readmatrix('local_fit_with_cebpa_runx1_general_PU1CEfeedback_PrepC_constRunx1_2N_maxV1_123456789101112_010100.txt');
% [~,I] = min(previous_results(:,end));
% best_pars = previous_results(I,1:end-1);
% initial_points = best_pars([1:17,23:26]);
% initial_points([4,18,19,21,22,23]) = initial_points([4,18,19,21,22,23])-0.01;
% initial_points([10,16,20]) = initial_points([10,16,20])+0.01;

lb = zeros(1,nvars)-6;
ub = zeros(1,nvars)+6;

% hill pars for feedback (Amin,Amax,N,K) % if using fitted HIll pars from
% no feedback, then '_constrainByFittedHill'
% lb(13:16) = best_pars(16:19);
% ub(13:16) = best_pars(16:19);
lb(21:24) = [0.0001,1,2,0.001];
ub(21:24) = [0.1,1,2,5];

lb(1:7) = -6;
ub(1:7) = 6;

% runx1_bcell_conc = p(18);
% runx1_ProLow_conc = p(19);
% runx1_ProHigh_conc = p(20);
lb(18:20) = [0,0,0];
ub(18:20) = [1,1,1];

% rescaling factor for promoter-only construct
% promoter_only_rescale_factor = p(17);
lb(17) = 0.001;
ub(17) = 5;

% nuc destab
lb(8:10) = 0;
ub(8:10) = 10;

% take into account whehteehr the energy exists in the first place. If not,
% upper bound is also 0
ub(8) = ub(8)*(nuc_destab_location(1)==1 || nuc_destab_location(2)==1);
ub(9) = ub(9)*(nuc_destab_location(3)==1 || nuc_destab_location(4)==1);
ub(10) = ub(10)*(nuc_destab_location(5)==1 || nuc_destab_location(6)==1);

% recruitmnet of RNAP
lb(11:13) = -6;
ub(11:13) = 0;

% runx pu14, runx cebp attraction 
lb(14:15) = -6;
ub(14:15) = 0;

% pu1 cebpa repulsion (set this to negative because I found evidence that
% they work together to open chromatin)
lb(16) = 0;
ub(16) = 10;

% particleswarm options
% options = optimoptions(@particleswarm,'UseParallel',false,'Display','iter','InitialSwarmMatrix',initial_points);
% %     'InitialSwarmMatrix',init_point);

% ga options
options = optimoptions('ga','Display','iter');

for itr = 1:10
% [x,fval] = particleswarm(fun,nvars,lb,ub);
% [x,fval] = particleswarm(fun,nvars,lb,ub,options);
[x,fval] = ga(fun,nvars,[],[],[],[],lb,ub,[],options);
create_parameter_file(final_fn, x, fval, 'a+')
end

