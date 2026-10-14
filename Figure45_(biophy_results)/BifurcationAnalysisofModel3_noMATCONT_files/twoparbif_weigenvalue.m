function [x_ap2f,s_ap2f,f_ap2f,x_ap2b,s_ap2b,f_ap2b] = twoparbif_weigenvalue(odefile,updated_ap,x_ap1,range_ap1,pars,lpn,MyContSetOptions)
% I modified this so that only one additional parmaeter is analyzed on top
% of alpha. If there are multiple parameters that can be analyzed on top of
% alpha (or any single parameter), will create a loop *outside* of this
% function

% this file also only continues from a single LP. If there are multiple LP,
% do a loop *outside* of this file

% rangeap1 is the array of all evaluated single parameter in one-par
% analysis

% but this file can do forward and backward simultaneously (instead of
% setting b outside to be 0 or 1 respectively)
%% set options for LP continuation
opt=contset;
% opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxStepsize',MyContSetOptions.maxstepsize);
opt=contset(opt,'MaxStepsize',MyContSetOptions.minstepsize);
opt=contset(opt,'MaxNumPoints',MyContSetOptions.maxpt);
opt=contset(opt,'Eigenvalues',1);

%% update parameter array and find the state variables evaluated at limit point
ap1 = updated_ap(1);
pars(ap1) = range_ap1(lpn); % update parameter array such that first ap is at LP

num_x = length(x_ap1(:,1))-1; % number of state variables
xnew = x_ap1(1:num_x,lpn); % state variables evaluated at parameter 1 at which LP happens

[x0 v0]=init_LP_LP(odefile,xnew,pars',updated_ap); % store initial condition (pars need to be converted into a column for init_LPLP to work

%% continue LP forward
opt=contset(opt,'Backward',0);
[x_ap2f,v,s_ap2f,h,f_ap2f]=cont(@limitpoint,x0,v0,opt);
% range_ap12f = x_ap2(num_x+1:num_x+2,:);
% stabf = test_stability(f_ap2(:,2:end-1));

%% continue LP backward
opt=contset(opt,'Backward',1);
[x_ap2b,v,s_ap2b,h,f_ap2b]=cont(@limitpoint,x0,v0,opt);
% range_ap12b = x_ap2(num_x+1:num_x+2,:);
% stabb = test_stability(f_ap2(:,2:end-1));

% %% plot forward and backward bifucation diagram (ap1 on y axis, ap2 on x axis)
% % plot forward
% plot(range_ap12f(2,:), range_ap12f(1,:),'Color',PlotOptions.colorf,'LineWidth',PlotOptions.lw,'DisplayName','Forward');
% hold on
% plot(range_ap12b(2,:), range_ap12b(1,:),'Color',PlotOptions.colorb,'LineWidth',PlotOptions.lw,'DisplayName','Backward');
% xlabel(updated_apn{2})
% ylabel(updated_apn{1})
% title(strcat(updated_apn{1},{' '},'vs',{' '},updated_apn{1},{' '},'two-par'))
% hold on; axis square; box on;
% xlim(PlotOptions.xL);
% ylim(PlotOptions.yL);
% set(gca,'XScale',PlotOptions.xsc,'YScale',PlotOptions.ysc,...
%         'XTickLabelRotation', 0,...
%         'TickLength',PlotOptions.tckL,'FontSize',PlotOptions.Fs,...
%         'FontName','Arial');
% %         'XTick',PlotOptions.xtk,'YTick',PlotOptions.ytk,...,
% end