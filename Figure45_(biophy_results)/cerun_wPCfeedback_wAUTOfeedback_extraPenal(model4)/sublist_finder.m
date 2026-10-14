function valid_sublist = sublist_finder(possible_state_label_persite,site_index)

count = site_index;

total_num_sites = length(possible_state_label_persite(1,:));
if count < total_num_sites % as long as we have not reached the rightmost site
    sublist = sublist_finder(possible_state_label_persite,count+1); % all possible labels that the next site to the right can have
    num_config_last_site = length(sublist(:,1)); % number of labels that the next site to the right can have
    num_config_curr_site = nnz(possible_state_label_persite(:,count)); % number of labels that the current site can have
    num_config_curr_site_idx = find(possible_state_label_persite(:,count)); % index of the labels of the current site
    
    valid_sublist = zeros(num_config_curr_site*num_config_last_site,total_num_sites-count+1); % create an empty array so that each possible labels
    % at the current site can be appaneded to all configurations from the
    % next site to the right, respectively
    
    for n = 1:num_config_curr_site
        valid_sublist((num_config_last_site*(n-1)+1):(num_config_last_site*n),:) = [repelem(possible_state_label_persite(num_config_curr_site_idx(n),count),num_config_last_site)',sublist];
    end
    
else
    valid_sublist = possible_state_label_persite(find(possible_state_label_persite(:,count)),count);
end