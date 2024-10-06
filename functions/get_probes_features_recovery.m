%function to extract probed states, probe options, and observed choices
%from data file - split up by features

function [probes, options, num_faulty_trials] = get_probes_features_recovery(experimental_data, feature_names, probe_type)

    %find exp probe trials
    trial_idx = find(strcmp(experimental_data.display, probe_type));
    
    %remove faulty trials/responses
    faulty_trials = find(abs(diff(trial_idx) - mean(diff(trial_idx))) > 20);
    if ~isempty(faulty_trials)
       trial_idx(faulty_trials) = [];
    end
    
    num_faulty_trials = numel(faulty_trials);
    
    %find probed states
    probes_tmp = experimental_data.Image(trial_idx);
    
    %map to numbers
    for k = 1:size(probes_tmp,1)
       probes(k,1) = find(strcmp(probes_tmp{k}(1), feature_names));
       probes(k,2) = find(strcmp(probes_tmp{k}(3), feature_names));
    end
    
    %find probe options
    options_tmp(:,1) = experimental_data.probeImageLeft(trial_idx);
    options_tmp(:,2) = experimental_data.probeImageRight(trial_idx);
    
    %map to numbers
    for k = 1:size(options_tmp,1)
       %left option         
       options(k,1) = find(strcmp(options_tmp{k,1}(1), feature_names));
       options(k,2) = find(strcmp(options_tmp{k,1}(3), feature_names));
       
       %right option         
       options(k,3) = find(strcmp(options_tmp{k,2}(1), feature_names));
       options(k,4) = find(strcmp(options_tmp{k,2}(3), feature_names));
    end