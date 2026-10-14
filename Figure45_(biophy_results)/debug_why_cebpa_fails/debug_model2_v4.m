function debug_model2_v4(num_Hill,LOCALORNOT,model_origin)
% can fixed energy pars lead to oopotimized Hill pars that still lead to
% stable steady states?
selec = 1:11;
nuc_destab_location = [0,1,0,1];

if LOCALORNOT == 1
    iend = 20;
    ext = 'local';
    options = optimoptions(@particleswarm,'Display','iter','UseParallel',false);
else
    iend=1;
task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
rng(task_id,'twister')
ext = 'NOTS';
options = optimoptions(@particleswarm,'Display','iter','UseParallel',true);
end

if num_Hill == 1
    func = @(p) ODE_objective_0618_withBestG_solvefor1HillPars(selec,p,nuc_destab_location,model_origin);
    nvars = 4;
    lb = [0.0001,1,20,0.1];
    ub = [0.2,1,20,1];
    fn = strcat(ext,'fixedG_oneHill_',model_origin,'Initial','.txt');
elseif num_Hill == 2
    func = @(p) ODE_objective_0618_withBestG_solvefor2HillPars(selec,p,nuc_destab_location,model_origin);
    nvars = 8;
    lb = [0.0001,1,4,0.1,0.0001,1,4,0.1];
    ub = [0.2,1,4,1,0.2,1,4,1];
    fn = strcat(ext,'fixedG_twoHill_',model_origin,'Initial','.txt');
end



for i =1:iend
    [x,fval] = particleswarm(func,nvars,lb,ub,options);
    create_parameter_file(fn,x,fval,'a+')
end
end