function new_sub_idx_all = recode_sub_idx(data_file)

    new_sub(:,1) = 1:numel(unique(data_file(:,1)));
    new_sub(:,2) = unique(data_file(:,1), 'stable');
    new_sub_idx_all = [];
    for i = 1:size(new_sub,1)
        new_idx = find(data_file(:,1) == new_sub(i,2));
        new_sub_idx = zeros(numel(new_idx),1) + new_sub(i,1);
        new_sub_idx_all = [new_sub_idx_all; new_sub_idx];
    end
            

end