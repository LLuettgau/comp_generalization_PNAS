function correct_trials = check_choice_trials(preproc_outputs, model_type)

if model_type == 1


elseif model_type >= 2
    
    all_choice_trials.choice_trials_exp_prior = [preproc_outputs.exp_probes_features preproc_outputs.exp_options_features ...
                               preproc_outputs.exp_probes_correct preproc_outputs.exp_trial_type];

    all_choice_trials.choice_trials_inf_prior = [preproc_outputs.inf_probes_features preproc_outputs.inf_options_features ...
                               preproc_outputs.inf_probes_correct preproc_outputs.inf_trial_type];


    all_choice_trials.choice_trials_exp_trans = [preproc_outputs.trans_exp_probes_features preproc_outputs.trans_exp_options_features ...
                               preproc_outputs.trans_exp_probes_correct preproc_outputs.trans_exp_trial_type];

    all_choice_trials.choice_trials_inf_trans = [preproc_outputs.trans_inf_probes_features preproc_outputs.trans_inf_options_features ...
                               preproc_outputs.trans_inf_probes_correct preproc_outputs.trans_inf_trial_type];

    sched_prior = preproc_outputs.prior_sched_features;
    sched_trans = preproc_outputs.trans_sched_features;

    %Get TMs of learning phases to validate
    for k = 1:2

        if k == 1
            sched = sched_prior;
        else
            sched = sched_trans;
        end
            
        state_visit_counts = zeros(preproc_outputs.dims_SF);
            
        for i = 1:size(sched,1)-1
            for jj = 1:2
                current_state = sched(i,jj);
                next_state = sched(i+1,jj);

                %keep count of specific transition
                state_visit_counts(current_state,next_state) = state_visit_counts(current_state,next_state) + 1;
            end
        end

        %make TM
        if k == 1
            all_mat_values = sum(state_visit_counts > 0,2);
            if any(sum(state_visit_counts > 0,2) > 2)
                matrix_idx = find(all_mat_values == max(all_mat_values));
                row_values = state_visit_counts(matrix_idx, :);

                nonzero_indices = find(row_values > 0);
                if ~isempty(nonzero_indices)
                    min_nonzero_value = min(row_values(nonzero_indices));
                    row_values(row_values == min_nonzero_value) = 0;
                end
                
                state_visit_counts(matrix_idx,:) = row_values;
                TM_prior = state_visit_counts ./ max(state_visit_counts,1);

            else
                TM_prior = state_visit_counts ./ max(state_visit_counts,1);
            end
        else
            TM_trans = state_visit_counts ./ max(state_visit_counts,1);
        end
    end
    
     correct_trials.prior_exp = get_transitions(TM_prior, all_choice_trials.choice_trials_exp_prior, 1);
     correct_trials.prior_inf = get_transitions(TM_prior, all_choice_trials.choice_trials_inf_prior, 1);

     correct_trials.trans_exp = get_transitions(TM_trans, all_choice_trials.choice_trials_exp_trans, 1);
     correct_trials.trans_inf = get_transitions(TM_trans, all_choice_trials.choice_trials_inf_trans, 1);
     
    end
end




