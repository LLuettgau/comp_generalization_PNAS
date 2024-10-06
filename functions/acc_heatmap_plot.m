function [] = acc_heatmap_plot(all_data, results, condition, probe_type)

    Labels_heatmap_prior = {'Prior: 4-cycle probe',  'Prior: 4-line probe', ...
                            'Prior: 6-line probe', 'Prior: 6-cycle probe'};
    
    Labels_heatmap_trans = {'Transfer: 4-cycle probe', 'Transfer: 6-cycle probe'};

    subset_data = all_data(results(:,8) == condition,:);
    subset_data = flip(sortrows(subset_data,4));
    subset_data = subset_data(:, [2:3 5:6]);

    imagesc(subset_data)
    colormap('hot');
    colorbar;
    clim([0.3, 1]);
    xticks(1:size(subset_data, 2)); 
    if condition == 1
        xticklabels({Labels_heatmap_prior{[1 3]} Labels_heatmap_trans{[1 2]}}); 
    else
        xticklabels({Labels_heatmap_prior{[2 4]} Labels_heatmap_trans{[1 2]}}); 
    end
    pbaspect([1 15 1])

    box off
    set(gca, 'FontSize', 16); 

    if condition == 1 
        if probe_type == 1
            ylabel('4-Cycle Prior', 'FontSize', 20);
            title('Experience Probes', 'FontSize', 20);
        elseif probe_type == 2
            title('Inference Probes', 'FontSize', 20);
        end
    elseif condition == 2
        if probe_type == 1
            ylabel('6-Cycle Prior', 'FontSize', 20);
            title('Experience Probes', 'FontSize', 20);
         elseif probe_type == 2
            title('Inference Probes', 'FontSize', 20);
        end
    end




end