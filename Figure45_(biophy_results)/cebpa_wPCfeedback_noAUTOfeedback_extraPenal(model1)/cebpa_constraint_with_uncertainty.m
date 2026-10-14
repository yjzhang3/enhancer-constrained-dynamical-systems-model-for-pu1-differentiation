clear 
clc
% define the uncertainty in cebpa level in B cell, prolow, and prohigh

%% B vs. mac
% proB_over_mono_mean = 4/44.2; % friedman paper, proB over monocyte
% Bcell_over_mono = 1/4; % densely interconnected paper, Bcella2-a4 over MONO1
% Bcell_over_mono_bound = [proB_over_mono_mean,Bcell_over_mono];

Bcell_over_mono_1 = 2^(0.8479-7.4597); % from haemosphere
Bcell_over_mono_2 = 4/44.2; % from friedman paper, proB over monocyte
Bcell_over_mono_bound = [Bcell_over_mono_1,Bcell_over_mono_2];

%% pro vs. mac (from haemosphere)
pro_over_mono_1 = 2^(6.3871-7.4597);
pro_over_mono_2 = 1; % not gerater than 1, no matter if you're looking at celltype or celllines
pro_over_mono_bound = [pro_over_mono_1,pro_over_mono_2];

save('Cebpa_ct_constraint.mat',"Bcell_over_mono_bound","pro_over_mono_bound")

%% pro vs. mac
% g1_mean = 100; % friedman paper, GMP
% g1_std = g1_mean-75.5;
% 
% m1_mean = 44.2; % friedman paper, monocyte
% m1_std = 51.6-m1_mean;
% 
% p1_low_mean = 0.34; %GATA1 paper, pu1low
% p1_low_std = 0.425-p1_low_mean;
% 
% p1_high_mean = 1.18; %GATA1 paper, pu1+
% p1_high_std = 1.28-p1_high_mean;
% 
% g2_mean = 2.1; % GATA1 paper, GMP
% g2_std = 2.21-g2_mean;
% 
% prolow_over_mono_mean = g1_mean/m1_mean*(p1_low_mean/g2_mean);
% prohigh_over_mono_mean = g1_mean/m1_mean*(p1_high_mean/g2_mean);
% 
% prolow_over_mono_std = error_propagation_fourMeasurements(g1_mean,m1_mean,p1_low_mean,g2_mean,g1_std,m1_std,p1_low_std,g2_std);
% prohigh_over_mono_std = error_propagation_fourMeasurements(g1_mean,m1_mean,p1_high_mean,g2_mean,g1_std,m1_std,p1_high_std,g2_std);
% 
% prolow_over_mono_bound = [prolow_over_mono_mean-prolow_over_mono_std,prolow_over_mono_mean+prolow_over_mono_std];
% prohigh_over_mono_bound = [prohigh_over_mono_mean-prohigh_over_mono_std,prohigh_over_mono_mean+prohigh_over_mono_std];

% save('Cebpa_ct_constraint.mat',"Bcell_over_mono_bound","prohigh_over_mono_bound","prolow_over_mono_bound")

% inline function
function X_std = error_propagation_fourMeasurements(g1_mean,m1_mean,p1_mean,g2_mean,g1_std,m1_std,p1_std,g2_std)

%X = g1/m1*(p1/g2);
delXdelg1 = p1_mean/(m1_mean*g2_mean);
delXdelm1 = g1_mean*p1_mean/g2_mean*(-m1_mean^(-2));
delXdelp1 = g1_mean/(m1_mean*g2_mean);
delXdelg2 = g1_mean*p1_mean/m1_mean*(-g2_mean^(-2));

X_std = sqrt(delXdelg1^2*g1_std^2 + delXdelm1^2*m1_std^2 + delXdelp1^2*p1_std^2 + delXdelg2^2*g2_std^2);
end