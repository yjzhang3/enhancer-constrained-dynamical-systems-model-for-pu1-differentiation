function fit_with_cebpa_runx1_general_PU1CEfeedback(objFuncHandle,nuc_destab_location,selec)

% selec_DD234567_F234_E23_AML = 1:12;
% selec = selec_DD234567_F234_E23_AML;
changePU1conc = 0;
load('error_prop_0216.mat')

fun = @(p) objFuncHandle(selec,p,changePU1conc,0,0,nuc_destab_location);
funname = func2str(objFuncHandle);

orig_name = strcat('local_fit_',funname(11:end),'_PattrC','_2N','_maxV1'); 

fn = strcat(orig_name,'_',num2str(selec),'_',num2str(nuc_destab_location));
fn = strrep(fn, ' ', '');
final_fn = strcat(fn,'.txt');
% task_fn = fullfile(NOTS_path,final_fn);
nvars = 26;

%% data constraint on concentration
pu1_bcell_conc_range = conc_ratio_all(1,[1,3]); 
pu1_ProLow_conc_range = conc_ratio_all(2,[1,3]); 
pu1_ProHigh_conc_range = conc_ratio_all(3,[1,3]); 

% cebp_bcell_conc_range = [4/43,4/43];
% cebp_ProHigh_conc_range = [100/43*1.18/2.1,100/43*1.18/2.1];
% cebp_ProLow_conc_range = [100/43*0.33/2.1,100/43*0.33/2.1];

% %%%%% choose one of the three from the following to determine RUNX1
% participation
% 1. distinct RUNX1 (most participation of runx1: distinct RUNX1 level across cell types, and
% runx1 participates in DNA binding, TF interaction, and TF-RNAP
% interaction)
% runx1_bcell_conc_range = [0.01,1];
% runx1_pro_conc_range = [0.18/0.22,1];

% 2. identical RUNX1 (less participation from runx1: identical runx2 level
% across cell types, and runx1 participates in DNA binding, TF interaction, and TF-RNAP
% interaction)
runx1_bcell_conc_range = [0,0.5];
runx1_pro_conc_range = [0.5,2];
% CD31+ly06C+(mac) has a median around 10^2 (let's say 100-500); Lin-c-kit+ has a median around
% 10^3 (let's say 300-700); B220+Cd43- has a median around 10^1 (let's say
% 8-10)


%% parameter bound
lb = zeros(1,nvars)-6;
ub = zeros(1,nvars)+6;
lb(1:7) = -6;
ub(1:7) = 6;

% p_cebp_hill = p(23:26); Amin,Amax,N,K
% lb(23:26) = [0,cebp_ProHigh_conc_range(1),1,0];
lb(23:26) = [0.0001,1,2,0.4];
ub(23:26) = [0.1,1,2,5];

% pu1_bcell_conc = p(18);
% pu1_ProLow_conc = p(19);
% pu1_ProHigh_conc = p(20);
lb(18:20) = [pu1_bcell_conc_range(1),pu1_ProLow_conc_range(1),pu1_ProHigh_conc_range(1)];
ub(18:20) = [pu1_bcell_conc_range(2),pu1_ProLow_conc_range(2),pu1_ProHigh_conc_range(2)];

% runx1_bcell_conc = p(21);
% runx1_pro_conc = p(22);
lb(21:22) = [runx1_bcell_conc_range(1),runx1_pro_conc_range(1)];
ub(21:22) = [runx1_bcell_conc_range(2),runx1_pro_conc_range(2)];

% rescaling factor for promoter-only construct
% promoter_only_rescale_factor = p(17);
lb(17) = 0;
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
lb(14:15) = -6;
ub(14:15) = 0;

% pu1 cebpa repulsion (set this to negative because I found evidence that
% they work together to open chromatin)
lb(16) = -6;
ub(16) = 0;

% particleswarm options
options = optimoptions(@particleswarm,'UseParallel',true,'Display','iter');
%     'InitialSwarmMatrix',init_point);

for itr = 1:1
% [x,fval] = particleswarm(fun,nvars,lb,ub);
[x,fval] = particleswarm(fun,nvars,lb,ub,options);
create_parameter_file(final_fn, x, fval, 'a+')
end

