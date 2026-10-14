function debug_why_CEAUTO_necessary_v2(LOCALORNOT)
% can fixed energy pars (from model 4) lead to oopotimized Hill pars that still lead to
% stable steady states?
selec = 1:12;
nuc_destab_location = [0,1,0,1,0,0];

if LOCALORNOT == 1
    iend = 5;
    ext = 'local';
    options = optimoptions(@particleswarm,'Display','iter','UseParallel',false);
else
    iend=1;
task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
rng(task_id,'twister')
ext = 'NOTS';
options = optimoptions(@particleswarm,'Display','iter','UseParallel',true);
end

    func = @(p)  ODE_objective_cerun_0618_fixedG_find1HillPars(selec,p,nuc_destab_location);

    nvars = 4;
    lb = [0.0001,1,15,0.1];
    ub = [0.2,1,15,1];
    % fn = strcat(ext,'fixedG_1Hill.txt');
    fn = strcat(ext,'fixedG_1Hill_rescaleLS1.txt');
% end

% options = optimoptions(@particleswarm,'Display','iter','UseParallel',true);

for i =1:iend
    [x,fval] = particleswarm(func,nvars,lb,ub,options);
    create_parameter_file(fn,x,fval,'a+')
end
end