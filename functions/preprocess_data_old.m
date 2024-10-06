function [experimental_data, prior_sched, prior_sched_prim, prior_sched_sec, trans_sched, trans_sched_prim, trans_sched_sec, all_exp_stim_names, all_trans_stim_names, dims_SF, dims_SR, ...
          exp_probes, exp_options, obs_exp_choices, ...
          inf_probes, inf_options, obs_inf_choices, ...
          trans_exp_probes, trans_exp_options, trans_obs_exp_choices, ...
          trans_inf_probes, trans_inf_options, trans_obs_inf_choices, ...
          exp_probes_features, exp_options_features, ...
          inf_probes_features, inf_options_features, ...
          trans_exp_probes_features, trans_exp_options_features, ...
          trans_inf_probes_features, trans_inf_options_features] = ...
          preprocess_data(task_data, task_begin_idx, task_end_idx, subj, feature_names, model_type)

    %extract individual data
    experimental_data = task_data(task_begin_idx(subj):task_end_idx(subj),:);
    
    disp(subj)
    %% find prior learning schedule

    prior_sched_idx = find(strcmp(experimental_data.display, 'learningTrial'));
    prior_sched_tmp = experimental_data.Image(prior_sched_idx);

    curr_state = {};
    all_prior_sched_tmp = [];

    for j = 1:size(prior_sched_tmp,1)
        if isempty(curr_state)
            curr_state = prior_sched_tmp(j);
            all_prior_sched_tmp = [all_prior_sched_tmp; curr_state];
        end

        if strcmp(prior_sched_tmp(j),curr_state) == 0
            curr_state = {};
        end

    end

    %remove erroneous duplicate of trials from one subject's data set
    if experimental_data.ParticipantPrivateID(1) == 6577525
        all_prior_sched_tmp(277:279) = [];
    end

    % get names of prior learning phase stimuli
    all_exp_stim_names_idx1 = find(strcmp(experimental_data.display, 'probeTrial'));
    all_exp_stim_names_idx2 = find(strcmp(experimental_data.display, 'inferenceTrial'));
    all_exp_stim_names(1:12,1) = unique(experimental_data.Image(all_exp_stim_names_idx1));
    all_exp_stim_names(13:24,1) = unique(experimental_data.Image(all_exp_stim_names_idx2));

    if model_type == 1 || model_type == 3
        %map sched to numbers
        for k = 1:size(all_prior_sched_tmp,1)
            prior_sched(k,1) = find(strcmp(all_prior_sched_tmp(k), all_exp_stim_names));
        end

        %create reshaped sched
        prior_sched = reshape(prior_sched,12,[])';
        prior_sched_tmp_reshaped = reshape(all_prior_sched_tmp,12,[])';

    end
    

    if model_type >= 2
        %map sched to numbers
        for k = 1:size(all_prior_sched_tmp,1)
            prior_sched_features(k,1) = find(strcmp(all_prior_sched_tmp{k}(1), feature_names));
            prior_sched_features(k,2) = find(strcmp(all_prior_sched_tmp{k}(3), feature_names));
        end

        %create reshaped sched
        prior_sched_prim = reshape(prior_sched_features(:,1),12,[])';
        prior_sched_sec = reshape(prior_sched_features(:,2),12,[])';
    end

    %% find transfer learning schedule

    trans_sched_idx = find(strcmp(experimental_data.display, 'transferTrial'));
    trans_sched_tmp = experimental_data.Image(trans_sched_idx);

    curr_state = {};
    all_trans_sched_tmp = [];

    for j = 1:size(trans_sched_tmp,1)
        if isempty(curr_state)
            curr_state = trans_sched_tmp(j);
            all_trans_sched_tmp = [all_trans_sched_tmp; curr_state];
        end

        if strcmp(trans_sched_tmp(j),curr_state) == 0
            curr_state = {};
        end

    end

    % get names of transfer learning phase stimuli
    all_trans_stim_names_idx1 = find(strcmp(experimental_data.display, 'transferProbeTrial'));
    all_trans_stim_names_idx2 = find(strcmp(experimental_data.display, 'transferInferenceTrial'));
    all_trans_stim_names(1:12,1) = unique(experimental_data.Image(all_trans_stim_names_idx1));
    all_trans_stim_names(13:24,1) = unique(experimental_data.Image(all_trans_stim_names_idx2));

    if model_type == 1 || model_type == 3
        %map sched to numbers
        for k = 1:size(all_trans_sched_tmp,1)
            trans_sched(k,1) = find(strcmp(all_trans_sched_tmp(k), all_trans_stim_names));
        end
        %create reshaped sched
        trans_sched = reshape(trans_sched,12,[])';
        trans_sched_tmp_reshaped = reshape(all_trans_sched_tmp,12,[])';

    end

    if model_type >= 2
        %map  sched to numbers
        for k = 1:size(all_trans_sched_tmp,1)
            trans_sched_features(k,1) = find(strcmp(all_trans_sched_tmp{k}(1), feature_names));
            trans_sched_features(k,2) = find(strcmp(all_trans_sched_tmp{k}(3), feature_names));
        end

        %create reshaped sched
        trans_sched_prim = reshape(trans_sched_features(:,1),12,[])';
        trans_sched_sec = reshape(trans_sched_features(:,2),12,[])';

    end


    if model_type == 1 || model_type == 3
        %% get dimensions for SR matrix
        dims_SR = size(all_exp_stim_names,1);

        %% experience probes (prior learning)
        %find exp probe trials
        probe_type = 'probeTrial';
        all_stim_names = all_exp_stim_names;
        [exp_probes, exp_options, obs_exp_choices, ...
            exp_num_faulty_trials] = get_probes(experimental_data, all_stim_names, probe_type);

        %% inference probes (prior learning)
        %find inference probe trials
        probe_type = 'inferenceTrial';
        [inf_probes, inf_options, ...
            obs_inf_choices, inf_num_faulty_trials] = get_probes(experimental_data, all_stim_names, probe_type);

        %% experience probes (transfer learning)
        %find transfer exp probe trials
        probe_type = 'transferProbeTrial';
        all_stim_names = all_trans_stim_names;
        [trans_exp_probes, trans_exp_options, ...
            trans_obs_exp_choices, trans_exp_num_faulty_trials] = get_probes(experimental_data, all_stim_names, probe_type);

        %% inference probes (transfer learning)
        %find transfer inference probe trials
        probe_type = 'transferInferenceTrial';
        [trans_inf_probes, trans_inf_options, ...
            trans_obs_inf_choices, trans_inf_num_faulty_trials] = get_probes(experimental_data, all_stim_names, probe_type);

        %% number of faulty trials - compounds
        all_faulty_trials_comp(subj,1) = exp_num_faulty_trials + inf_num_faulty_trials + ...
            trans_exp_num_faulty_trials + trans_inf_num_faulty_trials;


    end

    if model_type >= 2
        %% get dimensions for SF matrix
        dims_SF = size(feature_names,1);

        %% experience probes (prior learning)
        %find exp probe trials
        probe_type = 'probeTrial';
        [exp_probes_features, exp_options_features, ...
            obs_exp_choices, exp_num_faulty_trials_features] = get_probes_features(experimental_data, feature_names, probe_type);

        %% inference probes (prior learning)
        %find inference probe trials
        probe_type = 'inferenceTrial';
        [inf_probes_features, inf_options_features, ...
            obs_inf_choices, inf_num_faulty_trials_features] = get_probes_features(experimental_data, feature_names, probe_type);

        %% experience probes (transfer learning)
        %find transfer exp probe trials
        probe_type = 'transferProbeTrial';
        [trans_exp_probes_features, trans_exp_options_features, ...
            trans_obs_exp_choices, trans_exp_num_faulty_trials_features] = get_probes_features(experimental_data, feature_names, probe_type);

        %% inference probes (transfer learning)
        %find transfer inference probe trials
        probe_type = 'transferInferenceTrial';
        [trans_inf_probes_features, trans_inf_options_features, ...
            trans_obs_inf_choices, trans_inf_num_faulty_trials_features] = get_probes_features(experimental_data, feature_names, probe_type);

        %% number of faulty trials - features
        all_faulty_trials_features(subj,1) = exp_num_faulty_trials_features + inf_num_faulty_trials_features + ...
            trans_exp_num_faulty_trials_features + trans_inf_num_faulty_trials_features;


    end

    if model_type == 2 || model_type >= 4
        exp_probes = [];
        exp_options = [];
        inf_probes = [];
        inf_options = [];
        trans_exp_probes = [];
        trans_exp_options = [];
        trans_inf_probes = [];
        trans_inf_options = [];

        prior_sched = [];
        trans_sched = [];

        dims_SR = [];

    elseif model_type == 1
        exp_probes_features = [];
        exp_options_features = [];
        inf_probes_features = [];
        inf_options_features = [];
        trans_exp_probes_features = [];
        trans_exp_options_features = [];
        trans_inf_probes_features = [];
        trans_inf_options_features = [];

        prior_sched_prim = [];
        prior_sched_sec = [];

        trans_sched_prim = [];
        trans_sched_sec = [];
        
        dims_SF = [];
    end


end