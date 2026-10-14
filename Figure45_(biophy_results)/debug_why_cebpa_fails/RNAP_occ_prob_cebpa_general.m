function [prob,rows_with_100,BF_list] = RNAP_occ_prob_cebpa_general(config_list,energy_table,concentration_table,nuc_destab_location)

BF_list = zeros(length(config_list(:,1)),1);
for nn = 1:length(config_list(:,1))
    BF_list(nn,1) = BF_per_state_with_cebpa_general(config_list(nn,:),energy_table,concentration_table,nuc_destab_location);
end

rows_with_100 = find(any(config_list == 100, 2));  % logical vector
prob = sum(BF_list(rows_with_100,1))/sum(BF_list(:,1));

end