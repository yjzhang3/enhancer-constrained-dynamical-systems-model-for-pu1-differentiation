function model3_bifur_on_Hill_pars(ct,ap_idx,apn)
% ct = cell type (1=bcell,2=prolow,3=prohigh)
close all;init;
odefile = @AllMATCONTFuncHandle; % this is the ODE file for MATCONT to read
disp(ct)
disp(apn)

%% some parameters that are constant throughout simulation
A = readmatrix("fixedG_1Hill.txt");
[~,I] = min(A(:,end));
p_cebp_hill = A(I,1:end-1);
p_cebp_hill_Amin = p_cebp_hill(1);
p_cebp_hill_Amax = p_cebp_hill(2);
p_cebp_hill_N = p_cebp_hill(3);
p_cebp_hill_K = p_cebp_hill(4);
% p_auto_hill_Amin = p_auto_hill(1);
% p_auto_hill_Amax = p_auto_hill(2);
% p_auto_hill_N = p_auto_hill(3);
% p_auto_hill_K = p_auto_hill(4);
pars = [p_cebp_hill_Amin,p_cebp_hill_Amax,p_cebp_hill_N,p_cebp_hill_K,ct];

%% Set plotting options
fign = 1;
hp=figure(fign);
hp.Position = [250,445,1346,574];
PlotOptions.Fs = 16; 
% load('color_pal.mat','cols')
PlotOptions.tckL = [0.015,0.015];
PlotOptions.lw = 3;
PlotOptions.xsc = 'linear'; 
PlotOptions.ysc = 'linear';

ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);

hold(ax1, 'on');
hold(ax2, 'on');
%% set continuation options
% apn = 'cebphill_{K}'; % active parameter name
% ap_idx = 4; % select alpha as parameter
MyContSetOptions.maxstepsize = 0.01; % maximum step size
MyContSetOptions.maxpt = 50; % maximum number of points
MyContSetOptions.mintepsize = 0.00001; % min step size

%% test different continuation direction, and different initial ap value
% ap_initial_par_val_arr = [0.8];
search_dir_arr = [0,1];
OutDirName = strcat('Prediction_','[x]_vs_',apn,sprintf('_ct=%d',ct));


count = 0;
% for ii = 1:length(ap_initial_par_val_arr)
    for kk = 1:length(search_dir_arr)
        count = count + 1;
        [rangealpha,all_err] = setup_for_onepar_EP_continuation(ct,odefile,MyContSetOptions,...
        search_dir_arr(kk),ap_idx,pars,PlotOptions,count,ax1,ax2,OutDirName);

        if rangealpha ~= 0
            save(strcat("CT=",num2str(ct),"_ErrVSPar_",apn,"_Plot",num2str(count),'.mat'),'rangealpha','all_err')
        end
    end
% end

%% finalize figures
xlabel(ax1, apn);
ylabel(ax1, '[pu]');

xlabel(ax2, apn);
ylabel(ax2, '[ce]');

sgtitle(strcat('steady state vs ',{' '},apn,',',{' '},sprintf('ct=%d',ct)));

set(ax1,'XScale',PlotOptions.xsc,'YScale',PlotOptions.ysc,...
            'XTickLabelRotation', 0,...
            'TickLength',PlotOptions.tckL,'FontSize',PlotOptions.Fs,...
            'FontName','Arial');
    %         'XTick',PlotOptions.xtk,'YTick',PlotOptions.ytk,...,

    set(ax2,'XScale',PlotOptions.xsc,'YScale',PlotOptions.ysc,...
            'XTickLabelRotation', 0,...
            'TickLength',PlotOptions.tckL,'FontSize',PlotOptions.Fs,...
            'FontName','Arial');
    %         'XTick',PlotOptions.xtk,'YTick',PlotOptions.ytk,...,


% legend('+a,-pu1,back','sing','+a,+pu1,back','sing','-a,-pu1,forw','sing','-a,+pu1,forw','sing')
saveas(hp,strcat('[x]_vs_',apn,sprintf('_ct=%d',ct),'.fig'))