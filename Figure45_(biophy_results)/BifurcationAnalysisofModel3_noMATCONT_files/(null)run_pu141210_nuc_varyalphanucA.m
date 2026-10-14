function run_pu141210_nuc_varyalphanucA(Mp)
% this is the two-par analysis starting from the LP found from the alpha
% [pu1] bifurcation diagram
clc;close all
model_name = 'modelpu141210_nuc';
odefile = @PU1autoregulationeqns;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% backward search, high alpha, highpu1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyContSetOptions_p1.b = 1; 
x0 = 1;
alpha = 1;

%% set continuation options for ap1
ap = 16; % select alpha as parameter
apn = 'alpha';
MyContSetOptions_p1.maxstepsize = 0.01; % maximum step size
MyContSetOptions_p1.maxpt = 1000; % maximum number of points
MyContSetOptions_p1.mintepsize = 0.005; % min step size

%% Set plotting options for ap1
PlotOptions_p1.Fs = 16; 
load('color_pal.mat','cols')
PlotOptions_p1.tckL = [0.015,0.015];
PlotOptions_p1.lw = 3;
% PlotOptions.xsc = 'log'; PlotOptions.ysc = 'log';
% PlotOptions.xL = [1e-2,1e0]; 
% PlotOptions.yL = [1e-3,1e1];
PlotOptions_p1.xsc = 'linear'; PlotOptions_p1.ysc = 'linear';
PlotOptions_p1.xL = [0,1]; 
PlotOptions_p1.yL = [0,1];

%% define parameters
Mc = 0;K_CEBP = 0; nucA = 1; % model1 has no C/EBPa and does not change pu1-nuc interaction
load('cell_cycle_pars_Kueh.mat','C_max','C_min','Kc','D','N') % load cell cycle parmaeters based on the parameters used by Kueh
load('X=rep_nuc_best_pars.mat','best_pars')
mbp = mean(best_pars(:,[1,2,3,4,5,6,7,8]));
G_PU14 = mbp(1);
G_PU12 = mbp(2);
G_PU10 = mbp(3);
G_nuc = mbp(4);
G_1nuc = mbp(5);
G_2nuc = mbp(6);
G_3nuc = mbp(7);
G_PU1412 = mbp(8);
pars = [G_PU14,G_PU12,G_PU10,G_PU1412,G_nuc,G_1nuc,G_2nuc,G_3nuc,Mp,Mc,C_max,C_min,Kc,D,N,alpha,nucA];

%% Find equilibrium for PU1
display('Find equilibrium for PU1')
odeeqnfun = @(t,y) pu141210_nuc_ode_eqns(t,y,pars);
[t,y] = ode45(odeeqnfun,[0 1000],x0);
yi = y(end,:);

%% equilibrium continuation when changing alpha, and plot
[xalpha,salpha,rangealpha,falpha] = oneparbif_weigenvalue(odefile,yi,pars,ap,MyContSetOptions_p1);

ivar = 1;
fign = 1;
hp=figure(fign);
hp.Position = [153,42,1361,1074];
subplot(1,4,1)
plotonepar_wptName(xalpha,salpha,rangealpha,falpha,ivar,PlotOptions_p1,cols(4,:));
xlabel(apn);
ylabel('[pu1]');
title(strcat('one-par',{' '},'[pu1] steady state vs ',{' '},apn));

%% find [pu1] at which limit point happens
lpn_all = []; % here I'm already obtaining the index out of all par values, not the index of the struct at which LP happens
for i = 1:length(salpha)
    if salpha(i).label == 'LP'
        lpn_all = [lpn_all, salpha(i).index];
    end
end

%% set limit point continuation for ap2
updated_ap = [16,17]; updated_apn = {'alpha','nucA'};
MyContSetOptions_p2.maxpt = 5000;
MyContSetOptions_p2.maxstepsize = 0.01;
MyContSetOptions_p2.minstepsize = 0.005;

%% set plot options for ap2
PlotOptions_p2.Fs = 16; 
PlotOptions_p2.tckL = [0.015,0.015];
PlotOptions_p2.lw = 3;
% PlotOptions.xsc = 'log'; PlotOptions.ysc = 'log';
% PlotOptions.xL = [1e-2,1e0]; 
% PlotOptions.yL = [1e-3,1e1];
PlotOptions_p2.xsc = 'linear'; PlotOptions_p2.ysc = 'linear';
PlotOptions_p2.xL = [37.27,60]; 
PlotOptions_p2.yL = [0,1];
PlotOptions_p2.colorf = [0.47 0.25 0.80];
PlotOptions_p2.colorb = [0.25 0.80 0.54];

%% LP continuation and plot
x = xalpha;
range_ap1 = rangealpha;
for n = 1:length(lpn_all)
    subplot(1+length(lpn_all),1,1+n)
    lpn = lpn_all(n);
    [x_ap2f,s_ap2f,f_ap2f,x_ap2b,s_ap2b,f_ap2b] = twoparbif_weigenvalue(odefile,updated_ap,x,range_ap1,pars,lpn,MyContSetOptions_p2);
    plottwopar_wptName(x_ap2f,s_ap2f,f_ap2f,x_ap2b,s_ap2b,f_ap2b,updated_apn,PlotOptions_p2)
    sgtitle(sprintf('Mp=%d',Mp))
    saveas(hp,strcat(model_name,'_',updated_apn{1},'_vs_',updated_apn{2},sprintf('_Mp=%d',Mp),'.jpg'))
end
