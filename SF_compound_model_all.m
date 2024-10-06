%Successor feature model learning associations between compound stimuli 

function [energy, outputs] = SF_compound_model_all(inputs, params_to_opt)
   
    fixed_params = inputs.fixed_params;    
    dims_SR = inputs.dims_SR;
    simulation = inputs.simulation;    
    gamma = inputs.gamma; %discounting factor gamma in SF update, less updating of more distal states (less likely to visit states further away)     
    multi_lr = inputs.multi_lr;  
    multi_tau = inputs.multi_tau;  

    %schedules
    prior_sched = inputs.prior_sched;
    trans_sched = inputs.trans_sched;

    %probes - prior
    exp_probes = inputs.exp_probes;
    exp_options = inputs.exp_options;
    obs_exp_choices = inputs.obs_exp_choices;

    inf_probes = inputs.inf_probes;
    inf_options = inputs.inf_options;
    obs_inf_choices = inputs.obs_inf_choices;

    %probes - transfer
    trans_exp_probes = inputs.trans_exp_probes;
    trans_exp_options = inputs.trans_exp_options;
    trans_obs_exp_choices = inputs.trans_obs_exp_choices;

    trans_inf_probes = inputs.trans_inf_probes;
    trans_inf_options = inputs.trans_inf_options;
    trans_obs_inf_choices = inputs.trans_obs_inf_choices;

    %correct answers
    exp_probes_correct = inputs.exp_probes_correct;
    inf_probes_correct = inputs.inf_probes_correct;
    trans_exp_probes_correct = inputs.trans_exp_probes_correct;
    trans_inf_probes_correct = inputs.trans_inf_probes_correct;

    correct_answer_prior = [exp_probes_correct inf_probes_correct];
    correct_answer_trans = [trans_exp_probes_correct trans_inf_probes_correct];
     
    %get/specify parameters
    alpha = params_to_opt(1:multi_lr); %learning rate
    tau = params_to_opt(multi_lr+1:multi_lr+multi_tau); %softmax stochasticity parameter

    all_subpch = [];
    
    for task = 1:2 %repeat for prior and transfer learning task
        clear succ_prob
        clear inf_succ_prob
        clear inf_succ_prob_diff
        clear succ_prob_diff
        
        clear sched
        clear all_exp_probes
        clear all_exp_options
        clear all_exp_obs_choice
        
        clear all_inf_probes
        clear all_inf_options
        clear all_inf_obs_choice
        
        if task == 1
            sched = prior_sched;
            all_exp_probes = exp_probes;
            all_exp_options = exp_options;
            all_exp_obs_choice = obs_exp_choices;
            
            all_inf_probes = inf_probes;
            all_inf_options = inf_options;
            all_inf_obs_choice = obs_inf_choices;
            
        else
            sched = trans_sched;
            all_exp_probes = trans_exp_probes;
            all_exp_options = trans_exp_options;
            all_exp_obs_choice = trans_obs_exp_choices;
            
            all_inf_probes = trans_inf_probes;
            all_inf_options = trans_inf_options;
            all_inf_obs_choice = trans_obs_inf_choices;
        end
        
        %initialize successor representation matrix
        M = zeros(dims_SR,dims_SR);
        if multi_lr == 1            
            alpha_curr = alpha;
        else
            alpha_curr = alpha(task);
        end

        %initialize successor representation matrix
        sched_vec = reshape(sched.',1,[])';
        
        % sched_vec = repmat(1:12, 1, 36)';
        %% loop over trials/episodes
        iTrial = 0;
        
        for t = 1:(size(sched_vec,1)-1)
            if t == 1
                M(:,:,1) = M + eye(dims_SR(1));
            end

            M(:,:,t+1) = M(:,:,t);

            onehot = zeros(1,dims_SR(1));
            onehot(sched_vec(t)) = 1;

            M(sched_vec(t),:,t+1)= M(sched_vec(t),:,t) + alpha_curr * (onehot + gamma * M(sched_vec(t+1),:,t) - M(sched_vec(t),:,t));
  
            %% Probe trials
            if mod(t,size(sched,2)) == 0 || t == size(sched_vec,1)-1 %probes occur after every 12th trial
                % but we take the representation learned until then, not
                % considering the next transition (except last trial, which
                % happens at very end)
            
                if t == size(sched_vec,1)-1 
                    M_tmp = M(:,:,end);
                else
                    M_tmp = M(:,:,end-1);
                end
                
                iTrial = iTrial + 1;

                % if iTrial == 36 && task == 2
                %     figure
                %     imagesc(M_tmp)
                %     colormap('hot')
                %     clim([0, 1.2]);
                %     yticks([]);
                %     yticklabels({});
                %     xticks([]);
                %     xticklabels({});
                %     pbaspect([2 2 1])
                %     set(gcf,'color','w')
                % 
                %     exportgraphics(gcf, 'model_representation_compounds.png', 'Resolution', 1200);
                % end
                %% experience probes
                %Softmax function arbitrating between alternative 1) next feature combination in the trajectory
                %and alternative 2) some other feature combination in the trajectory
                
                %left option's SR prob
                succ_prob(iTrial,1) = M_tmp(all_exp_probes(iTrial),all_exp_options(iTrial,1));
                
                %right option's SR prob
                succ_prob(iTrial,2) = M_tmp(all_exp_probes(iTrial),all_exp_options(iTrial,2));
    
                %% inference probes
                %left option's SR prob
                inf_succ_prob(iTrial,1) = M_tmp(all_inf_probes(iTrial),all_inf_options(iTrial,1));
                
                %right option's SR prob
                inf_succ_prob(iTrial,2) = M_tmp(all_inf_probes(iTrial,1),all_inf_options(iTrial,2));
            

                % if mod(iTrial,10) == 0
                %     if iTrial == 10
                %        figure;
                %        t = tiledlayout(1,4)
                %     end
                % 
                %     nexttile            
                %     imagesc(M_tmp)
                % end
            end
        end
        
        for choice_type = 1:2 %experience or inference probe
            if multi_tau == 1            
                tau_curr = tau;
            else
                tau_curr = tau(choice_type);
            end
            
            %% turn values into softmax choice probs
            if choice_type == 1
                all_succ_prob = succ_prob;
                subpch = zeros(size(succ_prob,1),1);
                obs_choice = all_exp_obs_choice;
            else
                all_succ_prob = inf_succ_prob;
                subpch = zeros(size(succ_prob,1),1);
                obs_choice = all_inf_obs_choice;
            end
            
            succ_prob_diff = all_succ_prob(:,1) - all_succ_prob(:,2); %next state occupancy difference between left/right option
            p = 1./(1 + exp(-succ_prob_diff/tau_curr)); %softmax rule
            pch = [p 1-p]; %choice probabilities for left/right option

            if choice_type == 1 && task == 1
                all_subpch = [];

                exp_choices_prior = get_choices(pch, subpch, obs_choice, ...
                                                correct_answer_prior(:,1), simulation);
                all_subpch = [all_subpch; exp_choices_prior.subpch];

            elseif choice_type == 2 && task == 1

                inf_choices_prior = get_choices(pch, subpch, obs_choice, ...
                                                correct_answer_prior(:,2), simulation);
                all_subpch = [all_subpch; inf_choices_prior.subpch];

            elseif choice_type == 1 && task == 2

                exp_choices_trans = get_choices(pch, subpch, obs_choice, ...
                                                correct_answer_trans(:,1), simulation);
                all_subpch = [all_subpch; exp_choices_trans.subpch];

            elseif choice_type == 2 && task == 2

                inf_choices_trans = get_choices(pch, subpch, obs_choice, ...
                                                correct_answer_trans(:,2), simulation);
                all_subpch = [all_subpch; inf_choices_trans.subpch];

            end
        end
    end

    if simulation == 0
        num_all_choices = size(all_subpch(~isnan(all_subpch)),1);
    else
        num_all_choices = [];
    end

    
    %% Compute log likelihood of experience and inference choices in prior and transfer learning
    if multi_lr == 1
        if multi_tau == 1

            if alpha < 0 || alpha > 1 || tau < 0 || tau > 100
                energy = 9999999;
            else
                energy = -nansum(log(all_subpch)); %calculate the negative log likelihood of parameters given the choices
            end

        else

            if alpha < 0 || alpha > 1 || tau(1) < 0 || tau(1) > 100 || tau(2) < 0 || tau(2) > 100
                energy = 9999999;
            else
                energy = -nansum(log(all_subpch)); %calculate the negative log likelihood of parameters given the choices
            end

        end

    elseif multi_lr == 2
        if multi_tau == 1

            if alpha(1) < 0 || alpha(1) > 1 || alpha(2) < 0 || alpha(2) > 1 || tau < 0 || tau > 100
                energy = 9999999;
            else
                energy = -nansum(log(all_subpch)); %calculate the negative log likelihood of parameters given the choices
            end

        else

            if alpha(1) < 0 || alpha(1) > 1 || alpha(2) < 0 || alpha(2) > 1 || tau(1) < 0 || tau(1) > 100 || tau(2) < 0 || tau(2) > 100
                energy = 9999999;
            else
                energy = -nansum(log(all_subpch)); %calculate the negative log likelihood of parameters given the choices
            end

        end

    end

    if simulation == 1
        energy = [];
        energy = [];
    end

    % Create the outputs structure
    outputs.prior_chprob_exp = exp_choices_prior.chprob;
    outputs.prior_chprob_inf = inf_choices_prior.chprob;
    outputs.trans_chprob_exp = exp_choices_trans.chprob;
    outputs.trans_chprob_inf = inf_choices_trans.chprob;

    outputs.prior_sim_res_corr_exp = exp_choices_prior.sim_responses_corr;
    outputs.prior_sim_res_corr_inf = inf_choices_prior.sim_responses_corr;
    outputs.trans_sim_res_corr_exp = exp_choices_trans.sim_responses_corr;
    outputs.trans_sim_res_corr_inf = inf_choices_trans.sim_responses_corr;

    outputs.prior_sim_choices_exp = exp_choices_prior.sim_choices;
    outputs.prior_sim_choices_inf = inf_choices_prior.sim_choices;
    outputs.trans_sim_choices_exp = exp_choices_trans.sim_choices;
    outputs.trans_sim_choices_inf = inf_choices_trans.sim_choices;

    outputs.num_all_choices = num_all_choices;