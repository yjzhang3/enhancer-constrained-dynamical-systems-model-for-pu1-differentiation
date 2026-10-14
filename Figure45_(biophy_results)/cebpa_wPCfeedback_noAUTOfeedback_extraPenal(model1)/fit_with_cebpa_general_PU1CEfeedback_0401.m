function fit_with_cebpa_general_PU1CEfeedback_0401(objFuncHandle,nuc_destab_location,selec)

changePU1conc = 0;

fun = @(p) objFuncHandle(selec,p,changePU1conc,0,0,nuc_destab_location);
funname = func2str(objFuncHandle);
load(fullfile(parentdir(pwd, 1), 'Error propagation and Normalization', 'error_prop_0216.mat'));

% orig_name = strcat('local_fit_',funname(11:end),'_PrepC','_10N','_maxV1','_constByFitHill','_geneDataPen'); 
orig_name = strcat('local_fit_',funname(11:end),'_PattrC','_freeN','_maxV1'); 

fn = strcat(orig_name,'_',num2str(selec),'_',num2str(nuc_destab_location));
fn = strrep(fn, ' ', '');
final_fn = strcat(fn,'.txt');
% task_fn = fullfile(NOTS_path,final_fn);
nvars = 19;

%% data constraint on concentration
pu1_bcell_conc_range = [conc_ratio_all(1,1),conc_ratio_all(1,3)]; 
pu1_ProLow_conc_range = [conc_ratio_all(2,1),conc_ratio_all(2,3)]; 
pu1_ProHigh_conc_range = [conc_ratio_all(3,1),conc_ratio_all(3,3)]; 
% max_fold_change_ProHighLow = pu1_ProHigh_conc_range(2)/pu1_ProLow_conc_range(1);
% cebp_bcell_conc_range = [4/43,4/43];
% cebp_ProHigh_conc_range = [100/43*1.18/2.1,100/43*1.18/2.1];
% cebp_ProLow_conc_range = [100/43*0.33/2.1,100/43*0.33/2.1];

%% parameter bound
lb = zeros(1,nvars)-6;
ub = zeros(1,nvars)+6;
lb(1:6) = -6;
ub(1:6) = 6;

% pu1_bcell_conc = p(13);
% pu1_ProLow_conc = p(14);
% pu1_ProHigh_conc = p(15);
lb(13:15) = [pu1_bcell_conc_range(1),pu1_ProLow_conc_range(1),pu1_ProHigh_conc_range(1)];
ub(13:15) = [pu1_bcell_conc_range(2),pu1_ProLow_conc_range(2),pu1_ProHigh_conc_range(2)];

% hill pars for feedback (Amin,Amax,N,K) % if using fitted HIll pars from
% no feedback, then '_constrainByFittedHill'
% lb(16:19) = [0.115,1,10,0.336];
% ub(16:19) = [0.115,1,10,0.336];
lb(16:19) = [0.0001,1,10,0.1];
ub(16:19) = [0.1,1,15,1];

% rescaling factor for promoter-only construct
% promoter_only_rescale_factor = p(17);
lb(12) = 0;
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

% particleswarm options
options = optimoptions(@particleswarm,'UseParallel',false,'Display','iter');
%     'InitialSwarmMatrix',init_point);

for itr = 1:20
% [x,fval] = particleswarm(fun,nvars,lb,ub);
[x,fval] = particleswarm(fun,nvars,lb,ub,options);
create_parameter_file(final_fn, x, fval, 'a+')
end

