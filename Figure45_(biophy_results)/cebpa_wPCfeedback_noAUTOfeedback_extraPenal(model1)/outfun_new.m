function stop = outfun_new(optimValues,state,PCn,CEPROD)
    stop = false;

    % Use a fixed filename, or parameterize it if needed
    %filename = strcat('latest_swarm_',getenv('SLURM_ARRAY_TASK_ID'),'.mat');
    filename = strcat('latest_swarm_',getenv('SLURM_ARRAY_TASK_ID'),'ODENOTSfit_PS_fsolve_',sprintf('PCn%.0f%.0f_',PCn(1),PCn(2)),CEPROD,'.mat');
    switch state
        case 'iter'
            latest_swarm = optimValues.swarm; %#ok<NASGU>
            save(filename, 'latest_swarm');

        case 'done'
            % Final write (ensures last state is saved)
            latest_swarm = optimValues.swarm; %#ok<NASGU>
            save(filename, 'latest_swarm');
    end
end