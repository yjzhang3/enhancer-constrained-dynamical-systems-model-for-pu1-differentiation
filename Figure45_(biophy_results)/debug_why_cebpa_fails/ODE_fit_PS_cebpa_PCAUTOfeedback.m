function ODE_fit_PS_cebpa_PCAUTOfeedback(moreTime,objFuncHandle,nuc_destab_location,selec,PCn,AUTOn,CEPROD)

task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
rng(task_id,'twister')

%% housekeeping
% load('wRUNX1_noODEsim_ce_conc_pars_0401.mat','x0')
% x0(2,1) = 0.1734;
fun = @(p) objFuncHandle(selec,p,nuc_destab_location,CEPROD); %v5 and before
% fun = @(p) objFuncHandle(selec,p,nuc_destab_location,x0); 
% funname = func2str(objFuncHandle);

% orig_name = strcat('ODEfit_',funname(15:end-17)); 

% fn = strcat(orig_name,'_',num2str(selec),'_',num2str(nuc_destab_location));
% fn = strrep(fn, ' ', '');
final_fn = strcat('ODENOTSfit_PS_fsolve_',sprintf('PCn%.0f%.0f_',PCn(1),PCn(2)),sprintf('AUTOn%.0f%.0f_',AUTOn(1),AUTOn(2)),CEPROD,'.txt');
% task_fn = fullfile(NOTS_path,final_fn);
nvars = 20;

%% decdie whether there was timeout
if moreTime == 1
    filename = strcat('latest_swarm_',getenv('SLURM_ARRAY_TASK_ID'),'ODENOTSfit_PS_fsolve_',sprintf('PCn%.0f%.0f_',PCn(1),PCn(2)),sprintf('AUTOn%.0f%.0f',AUTOn(1),AUTOn(2)),CEPROD,'.mat');
    load(filename);
	initial_points = latest_swarm;
	options = optimoptions(@particleswarm,'Display','iter','InitialSwarmMatrix',initial_points,'UseParallel',true,'OutputFcn',@(optimValues,state) outfun_new(optimValues, state, PCn, AUTOn,CEPROD));
else
	options = optimoptions(@particleswarm,'Display','iter','UseParallel',true,'OutputFcn',@(optimValues,state) outfun_new(optimValues, state, PCn, AUTOn,CEPROD));
end

%% parameter bound
lb = zeros(1,nvars)-6;
ub = zeros(1,nvars)+6;
lb(1:6) = -6;
ub(1:6) = 6;

% hill pars for feedback (Amin,Amax,N,K) % if using fitted HIll pars from
% no feedback, then '_constrainByFittedHill'
% lb(13:16) = best_pars(16:19);
% ub(13:16) = best_pars(16:19);
lb(13:16) = [0.0001,1,PCn(1),0.1];
ub(13:16) = [0.2,1,PCn(2),0.9];
lb(17:20) = [0.0001,1,AUTOn(1),0.1];
ub(17:20) = [0.2,1,AUTOn(2),0.9];

% rescaling factor for promoter-only construct
% promoter_only_rescale_factor = p(17);
lb(12) = 1;
ub(12) = 10;

% nuc destab
lb(7:8) = 0;
ub(7:8) = 6;

% take into account whehteehr the energy exists in the first place. If not,
% upper bound is also 0
ub(7) = ub(8)*(nuc_destab_location(1)==1 || nuc_destab_location(2)==1);
ub(8) = ub(9)*(nuc_destab_location(3)==1 || nuc_destab_location(4)==1);

% recruitmnet of RNAP
lb(9:10) = -6;
ub(9:10) = 0;

% pu1 cebpa interaction 
lb(11) = -6;
ub(11) = 0;


for itr = 1:1
[x,fval] = particleswarm(fun,nvars,lb,ub,options);
create_parameter_file(final_fn, x, fval, 'a+')
end

