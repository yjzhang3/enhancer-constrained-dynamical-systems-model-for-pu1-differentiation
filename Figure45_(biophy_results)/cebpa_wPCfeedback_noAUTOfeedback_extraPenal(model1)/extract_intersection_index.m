function concise_ind_arr = extract_intersection_index(ind_arr,thres)

concise_ind_arr = [ind_arr(1)];
curr_arr_ind = ind_arr(1);
for i = 2:length(ind_arr)
    if ind_arr(i)-curr_arr_ind > thres
        concise_ind_arr = [concise_ind_arr,ind_arr(i)];
        curr_arr_ind = ind_arr(i);
    end
end