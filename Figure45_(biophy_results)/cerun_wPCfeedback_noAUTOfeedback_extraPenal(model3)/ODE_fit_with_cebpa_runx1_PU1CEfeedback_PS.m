function ODE_fit_with_cebpa_runx1_PU1CEfeedback_PS(moreTime,objFuncHandle,nuc_destab_location,selec,PCn,CEPROD)

%% housekeeping
% load('error_prop_0216.mat'); 
%previous_results = readmatrix('NOTS_fit_with_cebpa_runx1_general_PU1CEfeedback_0401_PattrC_2N_maxV1_123456789101112_010101.txt');
% load('wRUNX1_noODEsim_ce_conc_pars_0401.mat','x0')
% x0(2,1) = 0.1734;
task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
rng(task_id,'twister')


fun = @(p) objFuncHandle(selec,p,nuc_destab_location,CEPROD);
% funname = func2str(objFuncHandle);

%orig_name = strcat('ODENOTSfit_PS_',funname(15:end),'_PattrC','_2N','_maxV1'); 

%fn = strcat(orig_name,'_',num2str(selec),'_',num2str(nuc_destab_location));
%fn = strrep(fn, ' ', '');
%final_fn = strcat(fn,'.txt');

final_fn = strcat('ODENOTSfit_PS_fsolve_extraP_',getenv('SLURM_ARRAY_TASK_ID'),sprintf('PCn%.0f%.0f_',PCn(1),PCn(2)),CEPROD,'.txt');
% task_fn = fullfile(NOTS_path,final_fn);
nvars = 23;

%% parameter bound
if moreTime == 1
    filename = strcat('latest_swarm_',getenv('SLURM_ARRAY_TASK_ID'),'ODENOTSfit_PS_fsolve_',sprintf('PCn%.0f%.0f_',PCn(1),PCn(2)),CEPROD,'.mat');
	load(filename)
	initial_points = latest_swarm;
	options = optimoptions(@particleswarm,'Display','iter','InitialSwarmMatrix',initial_points,'UseParallel',true,'OutputFcn',@(optimValues, state) outfun_new(optimValues, state, PCn, CEPROD) );
else
	options = optimoptions(@particleswarm,'Display','iter','UseParallel',true,'OutputFcn',@(optimValues, state) outfun_new(optimValues, state, PCn, CEPROD) );
end

options.HybridFcn = @fmincon;
lb = zeros(1,nvars)-6;
ub = zeros(1,nvars)+6;

% hill pars for feedback (Amin,Amax,N,K) % if using fitted HIll pars from
% no feedback, then '_constrainByFittedHill'
% lb(13:16) = best_pars(16:19);
% ub(13:16) = best_pars(16:19);
lb(20:23) = [0.0001,1,PCn(1),0.001];
ub(20:23) = [0.1,1,PCn(2),2];

lb(1:7) = -6;
ub(1:7) = 6;

% runx1_bcell_conc = p(18);
% runx1_ProLow_conc = p(19);
% runx1_ProHigh_conc = p(20);
load('ct_TFPROD_constraint.mat',"run_Bcell_over_mac_bound","run_pro_over_mac_bound")
lb(18:19) = [run_Bcell_over_mac_bound(1),run_pro_over_mac_bound(1)];
ub(18:19) = [run_Bcell_over_mac_bound(2),run_pro_over_mac_bound(2)];

% rescaling factor for promoter-only construct
% promoter_only_rescale_factor = p(17);
lb(17) = 1;
ub(17) = 10;

% nuc destab
lb(8:10) = 0;
ub(8:10) = 6;

% take into account whehteehr the energy exists in the first place. If not,
% upper bound is also 0
ub(8) = ub(8)*(nuc_destab_location(1)==1 || nuc_destab_location(2)==1);
ub(9) = ub(9)*(nuc_destab_location(3)==1 || nuc_destab_location(4)==1);
ub(10) = ub(10)*(nuc_destab_location(5)==1 || nuc_destab_location(6)==1);

% recruitmnet of RNAP
lb(11:13) = -6;
ub(11:13) = 0;

% runx pu14, runx cebp attraction 
lb(14:15) = 0;
ub(14:15) = 0;

% pu1 cebpa 
lb(16) = -6;
ub(16) = 0;

% particleswarm options
%options = optimoptions(@particleswarm,'Display','iter','InitialSwarmMatrix',initial_points,'UseParallel',false,'OutputFcn',@outfun_new);
%options = optimoptions(@particleswarm,'Display','iter','UseParallel',true,'OutputFcn',@outfun_new);

for itr = 1:1
% [x,fval] = particleswarm(fun,nvars,lb,ub);
[x,fval] = particleswarm(fun,nvars,lb,ub,options);
create_parameter_file(final_fn, x, fval, 'a+')
end

