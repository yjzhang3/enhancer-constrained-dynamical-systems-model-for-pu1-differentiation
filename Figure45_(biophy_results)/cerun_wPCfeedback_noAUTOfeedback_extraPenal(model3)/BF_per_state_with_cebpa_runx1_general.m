function BF = BF_per_state_with_cebpa_runx1_general(config,energy_table,concentration_table,nuc_destab_location)

curr_config = config;
% -100 = nucleosome bound at promoter
% -14, -12, -10 = nucleosome bound at TFBS -14, -12, -10
% -1 = anysite unbound
% 12 = PU.1 bound at -12
% 10 = PU.1 bound at -10
% 100 = RNAP bound

% now enumrate all possibilities for the -14 
% 2 = CEBPA bound alone
% 3 = RUNX1 bound  alone (assume all three RUNX1 sites are occupied as long
% as RUNX1 binds)
% 14 = PU.1 bound alone at -14kb
% 214 = PU.1 and C/EBPa bound together
% 314 = (PU.1 alone bound at -14) AND (RUNX1 bound  alone)
% 23 = (RUNX1 bound  alone) AND (C/EBPa bound  alone)
% 2314 = (PU.1 alone bound at -14) AND (C/EBPa bound  alone)

%% binding energy
NUC_TF_count = sum(curr_config == -14) + sum(curr_config == -12) + sum(curr_config == -10); % number of nucleosome bound at TFBS
NUC_12_count = sum(curr_config == -12);
NUC_PRO_count = sum(curr_config == -100); % number of nucleosome bound at promoter
PU1_14_count = sum(curr_config == 14);
PU1_12_count = sum(curr_config == 12);
PU1_10_count = sum(curr_config == 10);
CEBP_alone_count = sum(curr_config == 2); 
RUNX_alone_count = sum(curr_config == 3);
RUNX_CEBP_at14_count = sum(curr_config == 23); 
RUNX_PU1_at14_count = sum(curr_config == 314); 
PU1_CEBP_at14_count = sum(curr_config == 214);
RUNX_PU1_CEBP_at14_count = sum(curr_config == 2314);
RNAP_count = sum(curr_config == 100);

% energy_table = [G_NUC,G_PU1_14,G_PU1_12,G_PU1_10,G_CEBP,G_RUNX,G_RNAP, ... %1-7
%     G_NUC_PU1,G_NUC_CEBP,G_NUC_RUNX1,... % 8-10
%     G_PU1_RNAP,G_CEBP_RNAP,G_RUNX1_RNAP,... % 11-13
%     G_RUNX_PU1,G_RUNX_CEBP,G_PU1_CEBP]; ... % 14-16

bind_e = (NUC_TF_count+NUC_PRO_count)*energy_table(1) + PU1_14_count*energy_table(2) + PU1_12_count*energy_table(3) + PU1_10_count*energy_table(4) + ...
    CEBP_alone_count*energy_table(5) + RUNX_alone_count*energy_table(6) + RNAP_count*energy_table(7) + ...
    RUNX_PU1_at14_count*(energy_table(2)+energy_table(6)) + RUNX_CEBP_at14_count*(energy_table(5)+energy_table(6)) + PU1_CEBP_at14_count*(energy_table(2)+energy_table(5)) + ... ...
    RUNX_PU1_CEBP_at14_count*(energy_table(2)+energy_table(5)+energy_table(6));

%% interaction energy
% nuc_destab_location =
% [PU1_PRO,PU1_allTFBS,C/EBPa_PRO,C/EBPa_12,RUNX1_PRO,RUNX1_12]

% energy_table = [G_NUC,G_PU1_14,G_PU1_12,G_PU1_10,G_CEBP,G_RUNX,G_RNAP, ... %1-7
%     G_NUC_PU1,G_NUC_CEBP,G_NUC_RUNX1,... % 8-10
%     G_PU1_RNAP,G_CEBP_RNAP,G_RUNX1_RNAP,... % 11-13
%     G_RUNX_PU1,G_RUNX_CEBP,G_PU1_CEBP]; ... % 14-16

PU1_PRO = nuc_destab_location(1);
PU1_allTFBS = nuc_destab_location(2);
CEBPa_PRO = nuc_destab_location(3);
CEBPa_12 = nuc_destab_location(4);
RUNX1_PRO = nuc_destab_location(5);
RUNX1_12 = nuc_destab_location(6);

total_PU1_14_count = PU1_14_count + RUNX_PU1_at14_count + RUNX_PU1_CEBP_at14_count + PU1_CEBP_at14_count;
total_CEBP_count = CEBP_alone_count + RUNX_CEBP_at14_count + PU1_CEBP_at14_count + RUNX_PU1_CEBP_at14_count;
total_RUNX_count = RUNX_alone_count + RUNX_PU1_at14_count + RUNX_CEBP_at14_count + RUNX_PU1_CEBP_at14_count;

int_e = PU1_PRO*(total_PU1_14_count+PU1_12_count+PU1_10_count)*(NUC_PRO_count)*energy_table(8)+...
    PU1_allTFBS*(total_PU1_14_count+PU1_12_count+PU1_10_count)*(NUC_TF_count)*energy_table(8)+...
    ...
    CEBPa_PRO*total_CEBP_count*NUC_PRO_count*energy_table(9)+...
    CEBPa_12*total_CEBP_count*NUC_12_count*energy_table(9)+...
    ...
    RUNX1_PRO*total_RUNX_count*NUC_PRO_count*energy_table(10)+...
    RUNX1_12*total_RUNX_count*NUC_12_count*energy_table(10)+...
    ...
    (total_PU1_14_count+PU1_12_count+PU1_10_count)*RNAP_count*energy_table(11)+...
    total_CEBP_count*RNAP_count*energy_table(12)+...
    total_RUNX_count*RNAP_count*energy_table(13)+...
    (RUNX_PU1_at14_count+RUNX_PU1_CEBP_at14_count)*energy_table(14)+... % RUNX1 attracts PU1 (only at -14kb)
    (RUNX_CEBP_at14_count+RUNX_PU1_CEBP_at14_count)*energy_table(15)+... % RUNX1 attracts CEBPa (only at -14kb)
    total_CEBP_count*(total_PU1_14_count+PU1_12_count+PU1_10_count)*energy_table(16);


curr_config_e = bind_e+int_e;

% concentration_table = [nuc],[pu1],[cebp],[runx],[RNAP]
curr_config_conc = concentration_table(1)^(NUC_TF_count + NUC_PRO_count) * concentration_table(2)^(total_PU1_14_count + PU1_12_count + PU1_10_count) ...
    * concentration_table(3)^(2*total_CEBP_count) ...
    * concentration_table(4)^(2*total_RUNX_count) * concentration_table(5)^RNAP_count;

%%
BF = curr_config_conc*exp(-curr_config_e);