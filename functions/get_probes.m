%function to extract probed states, probe options, and observed choices
%from data file

function [probes, options, correct_option, trial_type, obs_choices, num_faulty_trials] = get_probes(experimental_data, all_stim_names, probe_type)

    if ismember('Response', experimental_data.Properties.VariableNames)
        var_name = 'ScreenName';
    else
        var_name = 'display';
    end

    %find probe trials
    trial_idx = find(strcmp(experimental_data.(var_name), probe_type));
    
    %remove faulty trials/responses
    faulty_trials = find(abs(diff(trial_idx) - mean(diff(trial_idx))) > 20);
    if ~isempty(faulty_trials)
       trial_idx(faulty_trials) = [];
    end
    
    num_faulty_trials = numel(faulty_trials);
    
    %find participant choices
    if ismember('Response', experimental_data.Properties.VariableNames)
        obs_choices = experimental_data.Response(trial_idx);
    else
        obs_choices = [];
    end
    
    %find probed states
    probes_tmp = experimental_data.Image(trial_idx);
    
    %map to numbers
    for k = 1:size(probes_tmp,1)
       probes(k,1) = find(strcmp(probes_tmp(k), all_stim_names));
    end
    
    %find probe options
    options_tmp(:,1) = experimental_data.probeImageLeft(trial_idx);
    options_tmp(:,2) = experimental_data.probeImageRight(trial_idx);
    
    %map to numbers
    for k = 1:size(options_tmp,1)
       options(k,1) = find(strcmp(options_tmp(k,1), all_stim_names));
       options(k,2) = find(strcmp(options_tmp(k,2), all_stim_names));
    end

   %find correct option
   correct_option = experimental_data.correctAnswer(trial_idx);


   %find trial type 
   trial_type = experimental_data.featureProbe(trial_idx);