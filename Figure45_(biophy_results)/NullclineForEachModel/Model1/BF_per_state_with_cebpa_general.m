function BF = BF_per_state_with_cebpa_general(config,energy_table,concentration_table,nuc_destab_location)

curr_config = config;
% -100 = nucleosome bound at promoter
% -14, -12, -10 = nucleosome bound at TFBS -14, -12, -10
% -1 = anysite unbound
% 12 = PU.1 bound at -12
% 10 = PU.1 bound at -10
% 100 = RNAP bound

% now enumrate all possibilities for the -14 
% 2 = CEBPA bound alone
% as RUNX1 binds)
% 14 = PU.1 bound alone at -14kb
% 214 = PU.1 and C/EBPa bound together

%% binding energy
NUC_TF_count = sum(curr_config == -14) + sum(curr_config == -12) + sum(curr_config == -10); % number of nucleosome bound at TFBS
NUC_12_count = sum(curr_config == -12);
NUC_PRO_count = sum(curr_config == -100); % number of nucleosome bound at promoter
PU1_14_count = sum(curr_config == 14);
PU1_12_count = sum(curr_config == 12);
PU1_10_count = sum(curr_config == 10);
CEBP_alone_count = sum(curr_config == 2); 
PU1_CEBP_at14_count = sum(curr_config == 214);
RNAP_count = sum(curr_config == 100);

% energy_table = [G_NUC,G_PU1_14,G_PU1_12,G_PU1_10,G_CEBP,G_RNAP, ... %1-6
%     G_NUC_PU1,G_NUC_CEBP,... % 7-8
%     G_PU1_RNAP,G_CEBP_RNAP,... % 9-10
%     G_PU1_CEBP]; ... % 11

bind_e = (NUC_TF_count+NUC_PRO_count)*energy_table(1) + PU1_14_count*energy_table(2) + PU1_12_count*energy_table(3) + PU1_10_count*energy_table(4) + ...
    CEBP_alone_count*energy_table(5) +  RNAP_count*energy_table(6) + ...
    PU1_CEBP_at14_count*(energy_table(2)+energy_table(5));

%% interaction energy
% nuc_destab_location =
% [PU1_PRO,PU1_allTFBS,C/EBPa_PRO,C/EBPa_12,RUNX1_PRO,RUNX1_12]

% energy_table = [G_NUC,G_PU1_14,G_PU1_12,G_PU1_10,G_CEBP,G_RNAP, ... %1-6
%     G_NUC_PU1,G_NUC_CEBP,... % 7-8
%     G_PU1_RNAP,G_CEBP_RNAP,... % 9-10
%     G_PU1_CEBP]; ... % 11

PU1_PRO = nuc_destab_location(1);
PU1_allTFBS = nuc_destab_location(2);
CEBPa_PRO = nuc_destab_location(3);
CEBPa_12 = nuc_destab_location(4);

total_PU1_14_count = PU1_14_count + PU1_CEBP_at14_count;
total_CEBP_count = CEBP_alone_count + PU1_CEBP_at14_count;

int_e = PU1_PRO*(total_PU1_14_count+PU1_12_count+PU1_10_count)*(NUC_PRO_count)*energy_table(7)+...
    PU1_allTFBS*(total_PU1_14_count+PU1_12_count+PU1_10_count)*(NUC_TF_count)*energy_table(7)+...
    ...
    CEBPa_PRO*total_CEBP_count*NUC_PRO_count*energy_table(8)+...
    CEBPa_12*total_CEBP_count*NUC_12_count*energy_table(8)+...
    ...
    (total_PU1_14_count+PU1_12_count+PU1_10_count)*RNAP_count*energy_table(9)+...
    total_CEBP_count*RNAP_count*energy_table(10)+...
    total_CEBP_count*(total_PU1_14_count+PU1_12_count+PU1_10_count)*energy_table(11);


curr_config_e = bind_e+int_e;

% concentration_table = [nuc],[pu1],[cebp],[RNAP]
curr_config_conc = concentration_table(1)^(NUC_TF_count + NUC_PRO_count) * concentration_table(2)^(total_PU1_14_count + PU1_12_count + PU1_10_count) ...
    * concentration_table(3)^(2*total_CEBP_count) ...
    * concentration_table(4)^RNAP_count;

%%
BF = curr_config_conc*exp(-curr_config_e);