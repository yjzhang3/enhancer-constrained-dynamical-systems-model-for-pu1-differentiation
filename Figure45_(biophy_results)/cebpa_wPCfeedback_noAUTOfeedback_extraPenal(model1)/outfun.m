function [state,options,optchanged] = outfun(options,state,flag)
optchanged = false; % This function does not stop the solver
persistent psw_position_file
persistent hist
optchanged = false;
switch flag
    case 'init'
        hist(:,:,1) = state.Population;
    case 'iter'
        ibest = state.Best(end);
        ibest = find(state.Score == ibest,1,'last');
        bestx = state.Population(ibest,:);
%         bestf = ODE_objective_2D_with_cebpa_PU1CEfeedback(bestx);
        bestf = 0;
%         nextline = [optimValues.iteration,optimValues.bestx,optimValues.bestfval];
%         nextline = [state.Generation,bestx,bestf];
%         disp(state.iteration)
        hist = cat(3, hist, state.Population);
        psw_position_file = 'psw_position.txt';
        create_parameter_file(psw_position_file, bestx,bestf,'a+');
%         create_parameter_file(psw_position_file, optimValues.bestx, optimValues.bestfval,'a+');
    case 'done'
        assignin('base','hist',hist);
        assignin('base','psw_position_file',psw_position_file);
end
end