function all_construct_pu1_dependent_behavior(RNAPFuncHandle,energy_table,all_sublist_str,cebpa_mac,cebpa_bcell,runx1_mac,runx1_bcell,sim_d)
nuc_destab_location = [0,1,0,1,0,0];

pu1_totest_range = [0:0.01:10];

% sim_d = [D2_prob,D3_prob,D4_prob,D5_prob,D6_prob,D7_prob,F2_prob,F3_prob,F4_prob,E2_prob,E3_prob,AML1_prob];
D2_prob = sim_d(1);
D3_prob = sim_d(2);
D4_prob = sim_d(3);
D5_prob = sim_d(4);
D6_prob = sim_d(5);
D7_prob = sim_d(6);
F2_prob = sim_d(7);
F3_prob = sim_d(8);
F4_prob = sim_d(9);
E2_prob = sim_d(10);
E3_prob = sim_d(11);
AML1_prob = sim_d(12);

D2_valid_sublist = all_sublist_str.D2;
D3_valid_sublist = all_sublist_str.D3;
D4_valid_sublist = all_sublist_str.D4;
D5_valid_sublist = all_sublist_str.D5;
D6_valid_sublist = all_sublist_str.D6;
D7_valid_sublist = all_sublist_str.D7;
F2_valid_sublist = all_sublist_str.F2;
F3_valid_sublist = all_sublist_str.F3;
F4_valid_sublist = all_sublist_str.F4;
E2_valid_sublist = all_sublist_str.E2;
E3_valid_sublist = all_sublist_str.E3;
AML1_valid_sublist = all_sublist_str.AML1;

D2_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_mac,runx1_mac,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [D2_prob_unnorm, rows_with_100_D2, BF_list_D2] = RNAPFuncHandle(D2_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_D2, top2_BF_idx_D2] = maxk(BF_list_D2, 2);
            BF_sublist_D2 = BF_list_D2(rows_with_100_D2);
            [top2_BF_sub_D2, top2_BF_sub_idx_D2] = maxk(BF_sublist_D2, 2);
        else
            D2_prob_unnorm = RNAPFuncHandle(D2_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        D2_TR_range(i) = D2_prob_unnorm;
end

D3_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_mac,runx1_mac,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [D3_prob_unnorm, rows_with_100_D3, BF_list_D3] = RNAPFuncHandle(D3_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_D3, top2_BF_idx_D3] = maxk(BF_list_D3, 2);
            BF_sublist_D3 = BF_list_D3(rows_with_100_D3);
            [top2_BF_sub_D3, top2_BF_sub_idx_D3] = maxk(BF_sublist_D3, 2);
        else
            D3_prob_unnorm = RNAPFuncHandle(D3_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        D3_TR_range(i) = D3_prob_unnorm;
end

D4_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_mac,runx1_mac,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [D4_prob_unnorm, rows_with_100_D4, BF_list_D4] = RNAPFuncHandle(D4_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_D4, top2_BF_idx_D4] = maxk(BF_list_D4, 2);
            BF_sublist_D4 = BF_list_D4(rows_with_100_D4);
            [top2_BF_sub_D4, top2_BF_sub_idx_D4] = maxk(BF_sublist_D4, 2);
        else
            D4_prob_unnorm = RNAPFuncHandle(D4_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        D4_TR_range(i) = D4_prob_unnorm;
end

D5_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_mac,runx1_mac,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [D5_prob_unnorm, rows_with_100_D5, BF_list_D5] = RNAPFuncHandle(D5_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_D5, top2_BF_idx_D5] = maxk(BF_list_D5, 2);
            BF_sublist_D5 = BF_list_D5(rows_with_100_D5);
            [top2_BF_sub_D5, top2_BF_sub_idx_D5] = maxk(BF_sublist_D5, 2);
        else
            D5_prob_unnorm = RNAPFuncHandle(D5_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        D5_TR_range(i) = D5_prob_unnorm;
end

D6_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_mac,runx1_mac,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [D6_prob_unnorm, rows_with_100_D6, BF_list_D6] = RNAPFuncHandle(D6_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_D6, top2_BF_idx_D6] = maxk(BF_list_D6, 2);
            BF_sublist_D6 = BF_list_D6(rows_with_100_D6);
            [top2_BF_sub_D6, top2_BF_sub_idx_D6] = maxk(BF_sublist_D6, 2);
        else
            D6_prob_unnorm = RNAPFuncHandle(D6_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        D6_TR_range(i) = D6_prob_unnorm;
end

D7_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_mac,runx1_mac,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [D7_prob_unnorm, rows_with_100_D7, BF_list_D7] = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_D7, top2_BF_idx_D7] = maxk(BF_list_D7, 3);
            BF_sublist_D7 = BF_list_D7(rows_with_100_D7);
            [top2_BF_sub_D7, top2_BF_sub_idx_D7] = maxk(BF_sublist_D7, 3);

            disp('dominating configurations out of all configurations')
            disp(D7_valid_sublist(top2_BF_idx_D7,:)) % dominating configurations out of all configurations
            disp(top2_BF_D7)
            
            disp('dominating confgiurations out of all RNAP bound configurations')
            RNAPbound = D7_valid_sublist(rows_with_100_D7,:);
            disp(RNAPbound(top2_BF_sub_idx_D7,:)) % dominating confgiurations out of all RNAP bound configurations
            disp(top2_BF_sub_D7)
        else
            D7_prob_unnorm = RNAPFuncHandle(D7_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        D7_TR_range(i) = D7_prob_unnorm;
end

F2_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_mac,runx1_mac,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [F2_prob_unnorm, rows_with_100_F2, BF_list_F2] = RNAPFuncHandle(F2_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_F2, top2_BF_idx_F2] = maxk(BF_list_F2, 2);
            BF_sublist_F2 = BF_list_F2(rows_with_100_F2);
            [top2_BF_sub_F2, top2_BF_sub_idx_F2] = maxk(BF_sublist_F2, 2);
        else
            F2_prob_unnorm = RNAPFuncHandle(F2_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        F2_TR_range(i) = F2_prob_unnorm;
end

F3_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_mac,runx1_mac,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [F3_prob_unnorm, rows_with_100_F3, BF_list_F3] = RNAPFuncHandle(F3_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_F3, top2_BF_idx_F3] = maxk(BF_list_F3, 2);
            BF_sublist_F3 = BF_list_F3(rows_with_100_F3);
            [top2_BF_sub_F3, top2_BF_sub_idx_F3] = maxk(BF_sublist_F3, 2);
        else
            F3_prob_unnorm = RNAPFuncHandle(F3_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        F3_TR_range(i) = F3_prob_unnorm;
end

F4_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_mac,runx1_mac,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [F4_prob_unnorm, rows_with_100_F4, BF_list_F4] = RNAPFuncHandle(F4_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_F4, top2_BF_idx_F4] = maxk(BF_list_F4, 2);
            BF_sublist_F4 = BF_list_F4(rows_with_100_F4);
            [top2_BF_sub_F4, top2_BF_sub_idx_F4] = maxk(BF_sublist_F4, 2);
        else
            F4_prob_unnorm = RNAPFuncHandle(F4_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        F4_TR_range(i) = F4_prob_unnorm;
end

E2_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_bcell,runx1_bcell,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [E2_prob_unnorm, rows_with_100_E2, BF_list_E2] = RNAPFuncHandle(E2_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_E2, top2_BF_idx_E2] = maxk(BF_list_E2, 2);
            BF_sublist_E2 = BF_list_E2(rows_with_100_E2);
            [top2_BF_sub_E2, top2_BF_sub_idx_E2] = maxk(BF_sublist_E2, 2);
        else
            E2_prob_unnorm = RNAPFuncHandle(E2_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        E2_TR_range(i) = E2_prob_unnorm;
end

E3_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_bcell,runx1_bcell,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [E3_prob_unnorm, rows_with_100_E3, BF_list_E3] = RNAPFuncHandle(E3_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_E3, top2_BF_idx_E3] = maxk(BF_list_E3, 2);
            BF_sublist_E3 = BF_list_E3(rows_with_100_E3);
            [top2_BF_sub_E3, top2_BF_sub_idx_E3] = maxk(BF_sublist_E3, 2);
        else
            E3_prob_unnorm = RNAPFuncHandle(E3_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        E3_TR_range(i) = E3_prob_unnorm;
end

AML1_TR_range = zeros(size(pu1_totest_range));

for i = 1:length(pu1_totest_range)
        pu1 = pu1_totest_range(i);
        totest_concentration_table = [1,pu1,cebpa_mac,runx1_mac,1]; % [nuc],[pu1],[cebp],[runx],[RNAP]

        if i == length(pu1_totest_range)
            [AML1_prob_unnorm, rows_with_100_AML1, BF_list_AML1] = RNAPFuncHandle(AML1_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
            [top2_BF_AML1, top2_BF_idx_AML1] = maxk(BF_list_AML1, 2);
            BF_sublist_AML1 = BF_list_AML1(rows_with_100_AML1);
            [top2_BF_sub_AML1, top2_BF_sub_idx_AML1] = maxk(BF_sublist_AML1, 2);
        else
            AML1_prob_unnorm = RNAPFuncHandle(AML1_valid_sublist,energy_table,totest_concentration_table,nuc_destab_location);
        end

        AML1_TR_range(i) = AML1_prob_unnorm;
end

figure('position',[509,92,1094,906]); hold on;
colors = [
    0.00, 0.45, 0.70;  % blue
    0.85, 0.10, 0.10;  % red
    0.00, 0.65, 0.30;  % green
    0.90, 0.55, 0.00;  % orange
    0.55, 0.00, 0.75;  % purple
    0.00, 0.75, 0.80;  % cyan
    0.80, 0.40, 0.70;  % pink
    0.40, 0.25, 0.10;  % brown
    0.60, 0.75, 0.00;  % yellow-green
    0.95, 0.80, 0.10;  % yellow
    0.50, 0.50, 0.50;  % gray
    0.10, 0.20, 0.60;  % navy
];
constructs = {'D2','D3','D4','D5','D6','D7','F2','F3','F4','E2','E3','AML1'};
TR_ranges = {D2_TR_range,D3_TR_range,D4_TR_range,D5_TR_range,D6_TR_range,D7_TR_range, ...
             F2_TR_range,F3_TR_range,F4_TR_range,E2_TR_range,E3_TR_range,AML1_TR_range};

for k = 1:length(constructs)
    plot(pu1_totest_range, TR_ranges{k}, 'Color', colors(k,:), 'LineWidth', 3, 'DisplayName', constructs{k});
end

legend('show','location','southeast');
xlabel('PU.1 concentration');
ylabel('Normalized transcription rate');
set(gca,'FontSize',16,'Yscale','log')