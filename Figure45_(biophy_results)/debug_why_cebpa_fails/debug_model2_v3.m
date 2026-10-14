function debug_model2_v3(moreTime)
% when there is no feedback, what is the free cebpa for good pu1ss, pu1
% prod rate, and gene expression
% load("model2_best_results.mat",'best_pars','pu1_ss_conc','energy_table',"cebpa_ss_conc")
% promoter_only_rescale_factor = best_pars(12);
selec = 1:11;
nuc_destab_location = [0,1,0,1];

func = @(p) ODE_objective_0618_noPCfeedback(selec,p,nuc_destab_location);

task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
rng(task_id,'twister')

nvars = 15;
lb = zeros(1,nvars)-6;
ub = zeros(1,nvars)+6;
lb(1:6) = -6;
ub(1:6) = 6;

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

lb(13:15) = [0,0,0];
ub(13:15) = [1,1,1];

fn = 'noPCfeedback.txt';

%% decdie whether there was timeout
if moreTime == 1
    filename = strcat('latest_swarm_',getenv('SLURM_ARRAY_TASK_ID'),'.mat');
    load(filename);
	initial_points = latest_swarm;
	options = optimoptions(@particleswarm,'Display','iter','InitialSwarmMatrix',initial_points,'UseParallel',true,'OutputFcn',@outfun_new);
else
	options = optimoptions(@particleswarm,'Display','iter','UseParallel',true,'OutputFcn',@outfun_new);
end

for i =1:1
    [x,fval] = particleswarm(func,nvars,lb,ub,options);
    create_parameter_file(fn,x,fval,'a+')
end