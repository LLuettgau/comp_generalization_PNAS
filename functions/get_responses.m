%compute correct responses for different probe types

function [results, ...
          timed_out_trials, ...
          corr_resp_all_no_nan, ...
          corr_resp_all, ...
          corr_resp_prim, ...
          corr_resp_sec, ...
          all_corr_resp_prim, ...
          all_corr_resp_sec, ...
          rt_all, ...
          aggregated_data_file, ...
          aggregated_data_file_final, ...
          aggregated_data_file_first, ...
          aggregated_data_file_first_final, ...
          aggregated_data_file_both_phases, ...
          aggregated_data_file_both_phases_final, ...
          long_data_file_final] = get_responses(trial_num, ...
                                                experimental_data, ...
                                                trial_type, ...
                                                probe_type, ...
                                                results, ...
                                                timed_out_trials, ...
                                                corr_resp_all_no_nan, ...
                                                corr_resp_all, ...
                                                corr_resp_prim, ...
                                                corr_resp_sec, ...
                                                all_corr_resp_prim, ...
                                                all_corr_resp_sec, ...
                                                rt_all, ...
                                                aggregated_data_file_final, ...
                                                aggregated_data_file_first_final, ...
                                                aggregated_data_file_both_phases, ...
                                                aggregated_data_file_both_phases_final, ...
                                                long_data_file_final, ...
                                                curr_sub, ...
                                                col_idx)
                                                                  
                       
    %find probe trial idx
    trial_idx = find(strcmp(experimental_data.ScreenName, trial_type) ...
                             & strcmp(experimental_data.display, probe_type));
                          
                          
     if numel(trial_idx) ~= trial_num
       disp(['number probe trials: ' num2str(numel(trial_idx))])
    end
    
    %find time out trials
    timeout_idx = find(strcmp(experimental_data.ScreenName, trial_type) & experimental_data.TimedOut == 1);

    %corrected number of trials with a response
    trial_idx_corrected = trial_idx;
    trial_idx_nan = zeros(size(trial_idx,1),1);

    for j = 1:size(timeout_idx,1)
        trial_idx_corrected(trial_idx_corrected == timeout_idx(j)) = [];
        trial_idx_nan(trial_idx == timeout_idx(j)) = 1;
    end

    %percent correct during experience probe - excluding missing response trials 
    results(curr_sub,col_idx(1)) = nansum(experimental_data.Correct(trial_idx_corrected));
    results(curr_sub,col_idx(2)) = nansum(experimental_data.Correct(trial_idx_corrected)) / numel(trial_idx_corrected);
    

    %find different probe trial types during experience probe
    %primary probe trials
    primary_probes_trial_idx = find(strcmp(experimental_data.ScreenName, trial_type) & experimental_data.featureProbe == 1);
    primary_probes_trial_idx_corrected = primary_probes_trial_idx;
    
    %secondary probe trials
    secondary_probes_trial_idx = find(strcmp(experimental_data.ScreenName, trial_type) & experimental_data.featureProbe == 2);
    secondary_probes_trial_idx_corrected = secondary_probes_trial_idx;
    
    %number of trials with a response (corrected trials)
    primary_probes_trial_idx_nan = zeros(size(primary_probes_trial_idx,1),1);
    secondary_probes_trial_idx_nan = zeros(size(secondary_probes_trial_idx,1),1);

    for j = 1:size(timeout_idx,1)
        primary_probes_trial_idx_corrected(primary_probes_trial_idx_corrected == timeout_idx(j)) = [];
        primary_probes_trial_idx_nan(primary_probes_trial_idx == timeout_idx(j)) = 1;

        secondary_probes_trial_idx_corrected(secondary_probes_trial_idx_corrected == timeout_idx(j)) = [];
        secondary_probes_trial_idx_nan(secondary_probes_trial_idx == timeout_idx(j)) = 1;
    end
    
    
    %percent correct of different probe trial types during sequence
    %learning - excluding missing response trials 
    %primary
    results(curr_sub,col_idx(3)) = nansum(experimental_data.Correct(primary_probes_trial_idx_corrected));
    results(curr_sub,col_idx(4)) = nansum(experimental_data.Correct(primary_probes_trial_idx_corrected)) / numel(primary_probes_trial_idx_corrected);
    
    %secondary
    results(curr_sub,col_idx(5)) = nansum(experimental_data.Correct(secondary_probes_trial_idx_corrected));
    results(curr_sub,col_idx(6)) = nansum(experimental_data.Correct(secondary_probes_trial_idx_corrected)) / numel(secondary_probes_trial_idx_corrected);
    

    %check difference between total number and sum of NaNs in different
    %trial types in exp (should be 0)
    if strcmp(trial_type, 'probe') || strcmp(trial_type, 'transfer_probe')
        
        %write number of NaN trials to matrix
        timed_out_trials(curr_sub,2) = sum(trial_idx_nan);
        timed_out_trials(curr_sub,3) = (trial_num/2) - numel(primary_probes_trial_idx_corrected);
        timed_out_trials(curr_sub,4) = (trial_num/2) - numel(secondary_probes_trial_idx_corrected);
        timed_out_trials(curr_sub,6) = timed_out_trials(curr_sub,2) - sum(timed_out_trials(curr_sub,3:4));

    else
        %write number of NaN trials to matrix
        timed_out_trials(curr_sub,7) = sum(trial_idx_nan);
        timed_out_trials(curr_sub,9) = (trial_num/2) - numel(primary_probes_trial_idx_corrected);
        timed_out_trials(curr_sub,10) = (trial_num/2) - numel(secondary_probes_trial_idx_corrected);
        timed_out_trials(curr_sub,11) = timed_out_trials(curr_sub,7) - sum(timed_out_trials(curr_sub,9:10));
    end
    
    %write percent correct trajectories, overall and different trial types
    %- replace missing response trials with NaN
    %overall
    corr_resp_all_reduced = experimental_data.Correct(trial_idx)';
    corr_resp_all_no_nan(curr_sub,2:numel(trial_idx)+1) = corr_resp_all_reduced;
    corr_resp_all_reduced(trial_idx_nan==1) = NaN;

    corr_resp_all(curr_sub,2:numel(trial_idx)+1) = corr_resp_all_reduced;

    %primary
    corr_resp_prim_reduced = experimental_data.Correct(primary_probes_trial_idx)';
    corr_resp_prim_reduced(primary_probes_trial_idx_nan==1) = NaN;

    corr_resp_prim(curr_sub,2:numel(primary_probes_trial_idx)+1) = corr_resp_prim_reduced;
    
    %secondary
    corr_resp_sec_reduced = experimental_data.Correct(secondary_probes_trial_idx)';
    corr_resp_sec_reduced(secondary_probes_trial_idx_nan==1) = NaN;
    
    corr_resp_sec(curr_sub,2:numel(secondary_probes_trial_idx)+1) = corr_resp_sec_reduced;
    
    
    %write reaction time trajectories, overall
    %- replace missing response trials with NaN
    %overall
    rt_all_reduced = experimental_data.ReactionTime(trial_idx)' / 1000;
    rt_all_reduced(trial_idx_nan==1) = NaN;
    rt_all(curr_sub,2:numel(trial_idx)+1) = rt_all_reduced;
    
    
    %need primary/secondary timeseries with all trial entries and NaNs when trial
    %type was not presented on the current trial
    %primary
    all_trials_primary_probe_trials = NaN(1,size(trial_idx,1));

    all_trials_primary_probe_trials(experimental_data.featureProbe(trial_idx) == 1) = corr_resp_prim_reduced;
    all_corr_resp_prim(curr_sub,1) = experimental_data.ParticipantPrivateID(1);
    all_corr_resp_prim(curr_sub,2:numel(all_trials_primary_probe_trials)+1) = all_trials_primary_probe_trials;

    %secondary
    all_trials_secondary_probe_trials = NaN(1,size(trial_idx,1));
    
    all_trials_secondary_probe_trials(experimental_data.featureProbe(trial_idx) == 2) = corr_resp_sec_reduced;
    all_corr_resp_sec(curr_sub,1) = experimental_data.ParticipantPrivateID(1);
    all_corr_resp_sec(curr_sub,2:numel(all_trials_secondary_probe_trials)+1) = all_trials_secondary_probe_trials;

    %% aggregated data for first trials
    
    subidx = experimental_data.ParticipantPrivateID(1);
    
    aggregated_data_file_first(1:2,1) = subidx; %subID
    aggregated_data_file_first(1:2,2) = curr_sub; %subject number
    if sum(~isnan(all_corr_resp_prim(curr_sub,2:3))) == 0
        aggregated_data_file_first(1,3) = NaN; %primary correct
        aggregated_data_file_first(1,4) = NaN; %number of secondary trials
    else
        aggregated_data_file_first(1,3) = nansum(all_corr_resp_prim(curr_sub,2:3)); %primary correct
        aggregated_data_file_first(1,4) = sum(~isnan(all_corr_resp_prim(curr_sub,2:3))); %number of secondary trials
    end
    
    if sum(~isnan(all_corr_resp_sec(curr_sub,2:3))) == 0
        aggregated_data_file_first(2,3) = NaN; %secondary correct
        aggregated_data_file_first(2,4) = NaN; %number of secondary trials
    else
        aggregated_data_file_first(2,3) = nansum(all_corr_resp_sec(curr_sub,2:3)); %secondary correct
        aggregated_data_file_first(2,4) = sum(~isnan(all_corr_resp_sec(curr_sub,2:3))); %number of secondary trials
    end
    
    aggregated_data_file_first(1:2,5) = [-0.5 0.5];%trial types
    aggregated_data_file_first(1:2,6) = results(curr_sub,9);
    aggregated_data_file_first(1:2,7) = [1 2];%trial types
    aggregated_data_file_first(1:2,8) = results(curr_sub,8);
    
    %concatenate
    aggregated_data_file_first_final = [aggregated_data_file_first_final; aggregated_data_file_first];
    
    
    %% aggregated data for first trial comparison between prior and transfer - only use first occurences!
    end_trial = 2;
    
    aggregated_data_file_both_phases(1:2,1) = subidx; %subID
    aggregated_data_file_both_phases(1:2,2) = curr_sub; %subject number
    aggregated_data_file_both_phases(1:2,5) = [-0.5 0.5];%trial types
    aggregated_data_file_both_phases(1:2,6) = results(curr_sub,9);
    aggregated_data_file_both_phases(1:2,7) = [1 2];%trial types
    aggregated_data_file_both_phases(1:2,8) = results(curr_sub,8);
    

    if sum(~isnan(corr_resp_prim(curr_sub,2:end_trial))) == 0
        aggregated_data_file_both_phases(1,3) = NaN; %primary correct
        aggregated_data_file_both_phases(1,4) = NaN; %number of secondary trials
    else
        aggregated_data_file_both_phases(1,3) = nansum(corr_resp_prim(curr_sub,2:end_trial)); %primary correct
        aggregated_data_file_both_phases(1,4) = sum(~isnan(corr_resp_prim(curr_sub,2:end_trial))); %number of secondary trials
    end

    if sum(~isnan(corr_resp_sec(curr_sub,2:end_trial))) == 0
        aggregated_data_file_both_phases(2,3) = NaN; %secondary correct
        aggregated_data_file_both_phases(2,4) = NaN; %number of secondary trials
    else
        aggregated_data_file_both_phases(2,3) = nansum(corr_resp_sec(curr_sub,2:end_trial)); %secondary correct
        aggregated_data_file_both_phases(2,4) = sum(~isnan(corr_resp_sec(curr_sub,2:end_trial))); %number of secondary trials
    end

    %concatenate
    aggregated_data_file_both_phases_final = [aggregated_data_file_both_phases_final; aggregated_data_file_both_phases];

    
    
    
%     if strcmp(trial_type, 'probe') || strcmp(trial_type, 'probe_inf')
% 
%         if sum(~isnan(corr_resp_prim(curr_sub,2:end_trial))) == 0
%             aggregated_data_file_both_phases(1,3) = NaN; %primary correct
%             aggregated_data_file_both_phases(1,4) = NaN; %number of secondary trials
%         else
%             aggregated_data_file_both_phases(1,3) = nansum(corr_resp_prim(curr_sub,2:end_trial)); %primary correct
%             aggregated_data_file_both_phases(1,4) = sum(~isnan(corr_resp_prim(curr_sub,2:end_trial))); %number of secondary trials
%         end
% 
%         if sum(~isnan(corr_resp_sec(curr_sub,2:end_trial))) == 0
%             aggregated_data_file_both_phases(2,3) = NaN; %secondary correct
%             aggregated_data_file_both_phases(2,4) = NaN; %number of secondary trials
%         else
%             aggregated_data_file_both_phases(2,3) = nansum(corr_resp_sec(curr_sub,2:end_trial)); %secondary correct
%             aggregated_data_file_both_phases(2,4) = sum(~isnan(corr_resp_sec(curr_sub,2:end_trial))); %number of secondary trials
%         end
%         
%     elseif strcmp(trial_type, 'transfer_probe') || strcmp(trial_type, 'transfer_probe_inf')
%         if sum(~isnan(corr_resp_prim(curr_sub,2:end_trial))) == 0
%             aggregated_data_file_both_phases(1,5) = NaN; %primary correct
%             aggregated_data_file_both_phases(1,6) = NaN; %number of secondary trials
%         else
%             aggregated_data_file_both_phases(1,5) = nansum(corr_resp_prim(curr_sub,2:end_trial)); %primary correct
%             aggregated_data_file_both_phases(1,6) = sum(~isnan(corr_resp_prim(curr_sub,2:end_trial))); %number of secondary trials
%         end
% 
%         if sum(~isnan(corr_resp_sec(curr_sub,2:end_trial))) == 0
%             aggregated_data_file_both_phases(2,5) = NaN; %secondary correct
%             aggregated_data_file_both_phases(2,6) = NaN; %number of secondary trials
%         else
%             aggregated_data_file_both_phases(2,5) = nansum(corr_resp_sec(curr_sub,2:end_trial)); %secondary correct
%             aggregated_data_file_both_phases(2,6) = sum(~isnan(corr_resp_sec(curr_sub,2:end_trial))); %number of secondary trials
%         end
%         
%         %concatenate
%         aggregated_data_file_both_phases_final = [aggregated_data_file_both_phases_final; aggregated_data_file_both_phases];
%     end
    
    %% write long data file
    all_corr = experimental_data.Correct(trial_idx);
    
    %trial type indices
    trial_type_idx = experimental_data.featureProbe(trial_idx);
    
    for k = 1:size(trial_type_idx,1)
        if trial_type_idx == 1
            trial_type_idx_rec = -1;
        else
            trial_type_idx_rec = 1;
        end
    end   
        
    
    
    %aggregated data
    aggregated_data_file(1:2,1) = subidx; %subID
    aggregated_data_file(1:2,2) = curr_sub; %subject number
    aggregated_data_file(1,3) = results(curr_sub,col_idx(3)); %primary correct
    aggregated_data_file(1,4) = numel(primary_probes_trial_idx_corrected); %number of primary trials
    aggregated_data_file(2,3) = results(curr_sub,col_idx(5)); %secondary correct
    aggregated_data_file(2,4) = numel(secondary_probes_trial_idx_corrected); %number of secondary trials
    aggregated_data_file(1:2,5) = [-0.5 0.5];%trial types
    aggregated_data_file(1:2,6) = results(curr_sub,9);
    aggregated_data_file(1:2,7) = [1 2];%trial types
    aggregated_data_file(1:2,8) = results(curr_sub,8);
    
    %concatenate
    aggregated_data_file_final = [aggregated_data_file_final; aggregated_data_file];
    
    subidx = repmat(results(1),trial_num,1);
    %long data file
    long_data_file(1:trial_num,1) = subidx; %subID
    long_data_file(1:trial_num,2) = curr_sub; %subject number
    long_data_file(1:trial_num,3) = all_corr; %correct/incorrect response
    long_data_file(1:trial_num,4) = trial_type_idx_rec; %trial type
    long_data_file(:,5) = 1:size(long_data_file,1); %trial number
    long_data_file(:,6) = results(curr_sub,9);  %condition

    %concatenate
    long_data_file_final = [long_data_file_final; long_data_file];
                                                                                                                                              
end



