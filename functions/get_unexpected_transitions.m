%extract correct answers for specific transitions - 4-->1 (4-cycle) and 6-->1 (6-cycle) during
%transfer learning

function [agg_unexpected_transitions_tmp, first_unexpected_transitions_tmp, results] = get_unexpected_transitions(sched_raw, ...
                                                                             experimental_data, ...
                                                                             results, ...
                                                                             probes, ...
                                                                             probe_type_name, ...
                                                                             curr_sub, ...
                                                                             col_idx)

    
                                                                             
    %learning trials - Transfer
    all_learning_trial_trans_idx = find(strcmp(sched_raw.display, 'transferTrial'));
    all_learning_trials_trans = sched_raw.Image(all_learning_trial_trans_idx);
    
    %select those probe trials that were testing 4-->1 (4-cycle) and 6-->1
    %(6-cycle) transitions

    %define state names from trajectory
    %unexpected
    %4-->1
    four_cycle_41_transition_name = all_learning_trials_trans(4,1);
    %6-->1
    six_cycle_61_transition_name = all_learning_trials_trans(6,1);
    %4-->5 (if generalizing a 4 state cycle)
    six_cycle_45_transition_name = all_learning_trials_trans(5,1);
    
%     %expected transitions
%     %3-->4
%     four_cycle_34_transition_name = all_learning_trials_trans(3,1);
%     %5-->6
%     six_cycle_56_transition_name = all_learning_trials_trans(5,1);

    
    trans_trial_idx = find(strcmp(experimental_data.ScreenName, probes) ...
                             & strcmp(experimental_data.display, probe_type_name));
    
    all_probe_state_names = experimental_data.Image(trans_trial_idx);
    probe_type = experimental_data.featureProbe(trans_trial_idx);

    for k = 1:size(all_probe_state_names,1)
        all_probes(k,:) = [all_probe_state_names{k}(1) all_probe_state_names{k}(3) num2cell(trans_trial_idx(k)) probe_type(k)];
    end
    
    %find indices of probe trials probing those items (only in corresponding trial types!)
    four_cycle_41_idx = [];
    six_cycle_61_idx = [];
    six_cycle_45_idx = [];

    for k = 1:size(all_probes,1)
        if strcmp(four_cycle_41_transition_name{1}(1), all_probes{k,1}) && all_probes{k,4} == 1
           four_cycle_41_idx = [four_cycle_41_idx all_probes(k,3)];
        end
        
        if strcmp(six_cycle_61_transition_name{1}(3), all_probes{k,2}) && all_probes{k,4} == 2
            six_cycle_61_idx = [six_cycle_61_idx all_probes(k,3)];
        end
        
         if strcmp(six_cycle_45_transition_name{1}(3), all_probes{k,2}) && all_probes{k,4} == 2
            six_cycle_45_idx = [six_cycle_45_idx all_probes(k,3)];
        end
    end
    
%     %find specific expected 4 cycle and 6 cycle probes (for comparison)
%     expected_four_cycle_34_idx = [];
%     expected_six_cycle_56_idx = [];
% 
%     for k = 1:size(all_probes,1)
%         if strcmp(four_cycle_34_transition_name{1}(1), all_probes{k,1}) && all_probes{k,4} == 1
%            expected_four_cycle_34_idx = [expected_four_cycle_34_idx all_probes(k,3)];
%         end
%         
%         if strcmp(six_cycle_56_transition_name{1}(3), all_probes{k,2}) && all_probes{k,4} == 2
%             expected_six_cycle_56_idx = [expected_six_cycle_56_idx all_probes(k,3)];
%         end
%     end
%     
       
    %find all other 4 cycle and 6 cycle probes (for comparison)
    other_four_cycle_idx = [];
    other_six_cycle_idx = [];
    
    for k = 1:size(all_probes,1)
        if strcmp(four_cycle_41_transition_name{1}(1), all_probes{k,1}) == 0 && all_probes{k,4} == 1
           other_four_cycle_idx = [other_four_cycle_idx all_probes(k,3)];
        end
        
        if strcmp(six_cycle_61_transition_name{1}(3), all_probes{k,2}) == 0 && all_probes{k,4} == 2
            other_six_cycle_idx = [other_six_cycle_idx all_probes(k,3)];
        end
    end
    
    %check overlap between 4-cycle 4-->1 and all other 4-cycle transitions
    if numel(intersect(cell2mat(four_cycle_41_idx),cell2mat(other_four_cycle_idx))) > 0
        disp('overlap between four_cycle_41_idx and other_four_cycle_idx! Check that!')
    end
    
%     %check overlap between 4-cycle 4-->1 and expected 4-cycle transitions
%     if numel(intersect(cell2mat(four_cycle_41_idx),cell2mat(expected_four_cycle_34_idx))) > 0
%         disp('overlap between four_cycle_41_idx and other_four_cycle_idx! Check that!')
%     end


    %check overlap between 6-cycle 6-->1 and all other 6-cycle transitions
    if numel(intersect(cell2mat(six_cycle_61_idx),cell2mat(other_six_cycle_idx))) > 0
        disp('overlap between six_cycle_61_idx and other_six_cycle_idx! Check that!')
    end
     
%     %check overlap between 6-cycle 6-->1 and expected 6-cycle transitions
%     if numel(intersect(cell2mat(six_cycle_61_idx),cell2mat(expected_six_cycle_56_idx))) > 0
%         disp('overlap between expected_six_cycle_56_idx and six_cycle_61_idx! Check that!')
%     end
    
    subidx = experimental_data.ParticipantPrivateID(1);

    
    num_trials_to_include = 1;

    %compute percent correct for unexpected 4-->1 and 6-->1 transitions and expected transitions 
    if ~isempty(four_cycle_41_idx) && ~isempty(other_four_cycle_idx)
        four_cycle_41_correct = experimental_data.Correct(cell2mat(four_cycle_41_idx));
        results(curr_sub,col_idx(1)) = sum(four_cycle_41_correct) / size(four_cycle_41_correct,1);
        results(curr_sub,col_idx(2)) = mean(four_cycle_41_correct(1:num_trials_to_include));
        results(curr_sub,col_idx(17)) = experimental_data.trial_number(cell2mat(four_cycle_41_idx(1:num_trials_to_include)));
        
%         %compute percent correct for expected 4-cycle transitions
%         expected_four_cycle_34_correct = experimental_data.Correct(cell2mat(expected_four_cycle_34_idx));
%         results(curr_sub,col_idx(3)) = sum(expected_four_cycle_34_correct) / size(expected_four_cycle_34_correct,1);
%         results(curr_sub,col_idx(4)) = mean(expected_four_cycle_34_correct(1:num_trials_to_include));
         
        %compute percent correct for other 4-cycle
        other_four_cycle_correct = experimental_data.Correct(cell2mat(other_four_cycle_idx));
        results(curr_sub,col_idx(3)) = sum(other_four_cycle_correct) / size(other_four_cycle_correct,1);
        results(curr_sub,col_idx(4)) = mean(other_four_cycle_correct(1:num_trials_to_include));

        %compute overall accuracy for 4-cycle
        all_four_cycle_idx = [four_cycle_41_idx other_four_cycle_idx];
%         all_four_cycle_idx = [four_cycle_41_idx expected_four_cycle_34_idx];
        all_four_cycle_correct = experimental_data.Correct(cell2mat(all_four_cycle_idx));
        results(curr_sub,col_idx(5)) = sum(all_four_cycle_correct) / size(all_four_cycle_correct,1);

        
        
        %write to aggregate file
        agg_four_cycle_correct(1:2,1) = subidx; %subID
        agg_four_cycle_correct(1:2,2) = curr_sub; %subject number

        agg_four_cycle_correct(1,3) = sum(four_cycle_41_correct); %number of correct responses
        agg_four_cycle_correct(1,4) = size(four_cycle_41_correct,1); %number of 4-->1 trials
        agg_four_cycle_correct(1:2,5) = 1; %trial type
        agg_four_cycle_correct(1,6) = -0.5; %unexpected/expected

        agg_four_cycle_correct(2,3) = sum(other_four_cycle_correct); %number of correct responses
        agg_four_cycle_correct(2,4) = size(other_four_cycle_correct,1); %number of other four cycle probes
        agg_four_cycle_correct(2,6) = 0.5; %unexpected/expected
        
        agg_four_cycle_correct(1:2,7) = results(curr_sub,9); %condition
        agg_four_cycle_correct(1:2,8) = results(curr_sub,8); %condition recoded
        
        
        %first occurence of unexpected 4-->1 transition
        first_four_cycle_correct(1,1) = subidx; %subID
        first_four_cycle_correct(1,2) = curr_sub; %subject number

        first_four_cycle_correct(1,3) = nansum(four_cycle_41_correct(1:num_trials_to_include)); %number of correct responses
        first_four_cycle_correct(1,4) = sum(~isnan(four_cycle_41_correct(1:num_trials_to_include))); %number of 4-->1 trials
        first_four_cycle_correct(1,5) = 1; %trial type
        
        first_four_cycle_correct(1,6) = results(curr_sub,9); %condition
        first_four_cycle_correct(1,7) = results(curr_sub,8); %condition recoded
        
    else
        results(curr_sub,col_idx(1)) = NaN;
        results(curr_sub,col_idx(2)) = NaN;
        results(curr_sub,col_idx(3)) = NaN;
        results(curr_sub,col_idx(4)) = NaN;
        results(curr_sub,col_idx(5)) = NaN;
        results(curr_sub,col_idx(17)) = NaN;
        
        agg_four_cycle_correct = [];
        first_four_cycle_correct = [];
    end

   
    if ~isempty(six_cycle_61_idx) && ~isempty(other_six_cycle_idx)
        six_cycle_61_correct = experimental_data.Correct(cell2mat(six_cycle_61_idx));
        results(curr_sub,col_idx(6)) = sum(six_cycle_61_correct) / size(six_cycle_61_correct,1);
        results(curr_sub,col_idx(7)) = mean(six_cycle_61_correct(1:num_trials_to_include));
        results(curr_sub,col_idx(18)) = experimental_data.trial_number(cell2mat(six_cycle_61_idx(1:num_trials_to_include)));

%         %compute percent correct for other 6-cycle choices
%         expected_six_cycle_56_correct = experimental_data.Correct(cell2mat(expected_six_cycle_56_idx));
%         results(curr_sub,col_idx(8)) = sum(expected_six_cycle_56_correct) / size(expected_six_cycle_56_correct,1);
%         results(curr_sub,col_idx(9)) = mean(expected_six_cycle_56_correct(1:num_trials_to_include));
        
        %compute percent correct for other 6-cycle choices
        other_six_cycle_correct = experimental_data.Correct(cell2mat(other_six_cycle_idx));
        results(curr_sub,col_idx(8)) = sum(other_six_cycle_correct) / size(other_six_cycle_correct,1);
        results(curr_sub,col_idx(9)) = mean(other_six_cycle_correct(1:num_trials_to_include));
        
        %compute overall accuracy for 6-cycle
        all_six_cycle_idx = [six_cycle_61_idx other_six_cycle_idx];
%         all_six_cycle_idx = [six_cycle_61_idx expected_six_cycle_56_idx];
        all_six_cycle_correct = experimental_data.Correct(cell2mat(all_six_cycle_idx));
        results(curr_sub,col_idx(10)) = sum(all_six_cycle_correct) / size(all_six_cycle_correct,1);
        
        
        %write to aggregate file
        agg_six_cycle_61_correct(1:2,1) = subidx; %subID
        agg_six_cycle_61_correct(1:2,2) = curr_sub; %subject number

        agg_six_cycle_61_correct(1,3) = sum(six_cycle_61_correct); %number of correct responses
        agg_six_cycle_61_correct(1,4) = size(six_cycle_61_correct,1); %number of 6-->1 trials
        agg_six_cycle_61_correct(1:2,5) = 2; %trial type
        agg_six_cycle_61_correct(1,6) = -0.5; %unexpected/expected

        agg_six_cycle_61_correct(2,3) = sum(other_six_cycle_correct); %number of correct responses
        agg_six_cycle_61_correct(2,4) = size(other_six_cycle_correct,1); %number of other six cycle probes
        agg_six_cycle_61_correct(2,6) = 0.5; %unexpected/expected

        agg_six_cycle_61_correct(1:2,7) = results(curr_sub,9); %condition
        agg_six_cycle_61_correct(1:2,8) = results(curr_sub,8); %condition recoded

        %first occurence of unexpected 6-->1 transition
        first_six_cycle_61_correct(1,1) = subidx; %subID
        first_six_cycle_61_correct(1,2) = curr_sub; %subject number

        first_six_cycle_61_correct(1,3) = nansum(six_cycle_61_correct(1:num_trials_to_include)); %number of correct responses
        first_six_cycle_61_correct(1,4) = sum(~isnan(six_cycle_61_correct(1:num_trials_to_include))); %number of 6-->1 trials
        first_six_cycle_61_correct(1,5) = 2; %trial type
        
        first_six_cycle_61_correct(1,6) = results(curr_sub,9); %condition
        first_six_cycle_61_correct(1,7) = results(curr_sub,8); %condition recoded
        
    else
        results(curr_sub,col_idx(6)) = NaN;
        results(curr_sub,col_idx(7)) = NaN;
        results(curr_sub,col_idx(8)) = NaN;
        results(curr_sub,col_idx(9)) = NaN;
        results(curr_sub,col_idx(10)) = NaN;
        results(curr_sub,col_idx(18)) = NaN;

        agg_six_cycle_61_correct = [];
        first_six_cycle_61_correct = [];
    end
    
    
    
    
    %compute percent correct for 4-->5 choices
    if ~isempty(six_cycle_45_idx)
        six_cycle_45_correct = experimental_data.Correct(cell2mat(six_cycle_45_idx));
        results(curr_sub,col_idx(11)) = sum(six_cycle_45_correct) / size(six_cycle_45_correct,1);
        results(curr_sub,col_idx(12)) = mean(six_cycle_45_correct(1:num_trials_to_include));

        %write to aggregate file
        agg_six_cycle_45_correct(1,1) = subidx; %subID
        agg_six_cycle_45_correct(1,2) = curr_sub; %subject number

        agg_six_cycle_45_correct(1,3) = sum(six_cycle_45_correct); %number of correct responses
        agg_six_cycle_45_correct(1,4) = size(six_cycle_45_correct,1); %number of 4-->5 trials
        agg_six_cycle_45_correct(1,5) = 3; %trial type
        agg_six_cycle_45_correct(1,6) = -0.5; %unexpected

        agg_six_cycle_45_correct(1,7) = results(curr_sub,9); %condition
        agg_six_cycle_45_correct(1,8) = results(curr_sub,8); %condition recoded

        
        %first occurence of unexpected 4-->5 transition
        first_six_cycle_45_correct(1,1) = subidx; %subID
        first_six_cycle_45_correct(1,2) = curr_sub; %subject number

        first_six_cycle_45_correct(1,3) = nansum(six_cycle_45_correct(1:num_trials_to_include)); %number of correct responses
        first_six_cycle_45_correct(1,4) = sum(~isnan(six_cycle_45_correct(1:num_trials_to_include))); %number of 4-->5 trials
        first_six_cycle_45_correct(1,5) = 3; %trial type
        
        first_six_cycle_45_correct(1,6) = results(curr_sub,9); %condition
        first_six_cycle_45_correct(1,7) = results(curr_sub,8); %condition recoded
        
        
    else
        results(curr_sub,col_idx(11)) = NaN;
        results(curr_sub,col_idx(12)) = NaN;
        results(curr_sub,col_idx(13)) = NaN;
        results(curr_sub,col_idx(14)) = NaN;
        results(curr_sub,col_idx(15)) = NaN;
        
        agg_six_cycle_45_correct = [];
        first_six_cycle_45_correct = [];

    end
    

    %concatenate
    agg_unexpected_transitions_tmp = [agg_four_cycle_correct; agg_six_cycle_61_correct; agg_six_cycle_45_correct];
    first_unexpected_transitions_tmp = [first_four_cycle_correct; first_six_cycle_61_correct; first_six_cycle_45_correct];

end