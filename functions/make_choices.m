function [sim_res, sim_choices, w] = make_choices(M, SF_M, tau, w, omega, choice_trials, mechanism, trial_type, iTrial)
     
    if mechanism < 3 || mechanism == 4
        if trial_type == 1 %exp probes
            probe_state_opt = choice_trials.probes_exp;
            corr_opt = choice_trials.probes_exp_correct;
            incorr_opt = choice_trials.probes_exp_incorrect;
        else %inf probes
            probe_state_opt = choice_trials.probes_inf;
            corr_opt = choice_trials.probes_inf_correct;
            incorr_opt = choice_trials.probes_inf_incorrect;
        end
    else
        if trial_type == 1 %exp probes
            %compounds
            probe_state_SR_opt = choice_trials.SR.probes_exp;
            corr_SR_opt = choice_trials.SR.probes_exp_correct;
            incorr_SR_opt = choice_trials.SR.probes_exp_incorrect;
            
            %features
            probe_state_SF_opt = choice_trials.SF.probes_exp;
            corr_SF_opt = choice_trials.SF.probes_exp_correct;
            incorr_SF_opt = choice_trials.SF.probes_exp_incorrect;
            
        else %inf probes
            %compounds
            probe_state_SR_opt = choice_trials.SR.probes_inf;
            corr_SR_opt = choice_trials.SR.probes_inf_correct;
            incorr_SR_opt = choice_trials.SR.probes_inf_incorrect;
            
            %features
            probe_state_SF_opt = choice_trials.SF.probes_inf;
            corr_SF_opt = choice_trials.SF.probes_inf_correct;
            incorr_SF_opt = choice_trials.SF.probes_inf_incorrect;
            
        end

    end
    
    choice_opt_corr = [];
    choice_opt_incorr = [];
    
    %make choices
    if mechanism == 1 %SF compound model
  
        %correct option's SR prob
        succ_prob(1) = M(probe_state_opt(iTrial),corr_opt(iTrial));
        
        %incorrect option's SR prob
        succ_prob(2) = M(probe_state_opt(iTrial),incorr_opt(iTrial));
        
    elseif mechanism == 0 || mechanism == 2 || mechanism == 4 %SF feature model
        
        %calculate vector norm for presented options/successor features
        choice_opt_corr = [M(probe_state_opt(iTrial,1),corr_opt(iTrial,:)); ...
                           M(probe_state_opt(iTrial,2),corr_opt(iTrial,:))];
                       
        choice_opt_corr = sum(choice_opt_corr);
        
        succ_prob(1) = norm(choice_opt_corr);

        choice_opt_incorr = [M(probe_state_opt(iTrial,1),incorr_opt(iTrial,:)); ...
                             M(probe_state_opt(iTrial,2),incorr_opt(iTrial,:))];
                         
                         
        choice_opt_incorr = sum(choice_opt_incorr);

        succ_prob(2) = norm(choice_opt_incorr);
        
        if mechanism == 4 %weighting of representations
            succ_prob_trans = succ_prob;
            clear succ_prob

            %calculate vector norm for presented options/successor features
            choice_opt_corr_prior = [SF_M(probe_state_opt(iTrial,1),corr_opt(iTrial,:)); ...
                                     SF_M(probe_state_opt(iTrial,2),corr_opt(iTrial,:))];
                           
            choice_opt_corr_prior = sum(choice_opt_corr_prior);
            
            succ_prob_prior(1) = norm(choice_opt_corr_prior);
    
            choice_opt_incorr_prior = [SF_M(probe_state_opt(iTrial,1),incorr_opt(iTrial,:)); ...
                                       SF_M(probe_state_opt(iTrial,2),incorr_opt(iTrial,:))];
                             
                             
            choice_opt_incorr_prior = sum(choice_opt_incorr_prior);
    
            succ_prob_prior(2) = norm(choice_opt_incorr_prior);
            
            %weighting of prior learning and transfer learning SF representation
            %left option
            succ_prob(1) = w * succ_prob_prior(1) + (1-w) * succ_prob_trans(1);
            
            %right option
            succ_prob(2) = w * succ_prob_prior(2) + (1-w) * succ_prob_trans(2);
        end
        
    elseif mechanism == 3 %SF compound/SF feature hybrid 

        
        %correct option's SR prob
        succ_prob_SR(1) = M(probe_state_SR_opt(iTrial),corr_SR_opt(iTrial));
        
        %incorrect option's SR prob
        succ_prob_SR(2) = M(probe_state_SR_opt(iTrial),incorr_SR_opt(iTrial));
        
        
        %calculate vector norm for presented options/successor features
        %correct option
        choice_opt_corr = [SF_M(probe_state_SF_opt(iTrial,1),corr_SF_opt(iTrial,:)); ...
                           SF_M(probe_state_SF_opt(iTrial,2),corr_SF_opt(iTrial,:))];
                       
        choice_opt_corr = sum(choice_opt_corr);
        
        succ_prob_SF(1) = norm(choice_opt_corr);

        %incorrect option
        choice_opt_incorr = [SF_M(probe_state_SF_opt(iTrial,1),incorr_SF_opt(iTrial,:)); ...
                             SF_M(probe_state_SF_opt(iTrial,2),incorr_SF_opt(iTrial,:))];
                         
        choice_opt_incorr = sum(choice_opt_incorr);

        succ_prob_SF(2) = norm(choice_opt_incorr);
        
        
        %weighting of SR and SF representation
        %left option
        succ_prob(1) = w * succ_prob_SF(1) + (1-w) * succ_prob_SR(1);
        
        %right option
        succ_prob(2) = w * succ_prob_SF(2) + (1-w) * succ_prob_SR(2);

        if trial_type == 2
            succ_prob_diff_SR = succ_prob_SR(1) - succ_prob_SR(2);
            
            %softmax for SR probs
            inf_SRSF_p(1) = 1./(1 + exp(-succ_prob_diff_SR/tau));
            
            
            succ_prob_diff_SF = succ_prob_SF(1) - succ_prob_SF(2);
            
            %softmax for SF probs
            inf_SRSF_p(2) = 1./(1 + exp(-succ_prob_diff_SF/tau));
            
        end
       
    end
    
    
    %% Softmax
    %difference between successor state probabilities
    succ_prob_diff = succ_prob(1) - succ_prob(2); %next state occupancy difference between left/right option

    % turn successor state probabilities into softmax choice probs
    p = 1./(1 + exp(-succ_prob_diff/tau)); %softmax rule - probability of selecting left option
    pch = [p 1-p]; %choice probabilities for left/right option

    rand_choices = rand(1,1);

    sim_choices(rand_choices<pch(1)) = 1;
    sim_choices(rand_choices>pch(1)) = 2;

    chosen_option_correct = sim_choices == 1;
    sim_res = chosen_option_correct;

    %updating of weighting parameter
    if mechanism == 3 && trial_type == 2
        %update weighting parameter based on "confidence"/preference
        %strength in representations
        if inf_SRSF_p(1) >= inf_SRSF_p(2)
            w = w - omega * inf_SRSF_p(1);
        elseif inf_SRSF_p(1) < inf_SRSF_p(2)
            w = w + omega * inf_SRSF_p(2);
        end

        %prevent weighting parameter from growing infinitely
        if w > 1
            w = 1;
        elseif w < 0
            w = 0;
        end

    end

end