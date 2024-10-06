%Successor feature model learning associations between features with
%transfer learning (permutations)

function [energy, outputs] = SF_feature_model_all(inputs, params_to_opt)
             
    fixed_params = inputs.fixed_params;
    dims_SF = inputs.dims_SF;
    condition = inputs.condition;
    gamma = inputs.gamma; %discounting factor gamma in SF update, less updating of more distal states (less likely to visit states further away)     
    simulation = inputs.simulation;    
    transfer_model = inputs.transfer_model;  
    multi_lr = inputs.multi_lr;  
    multi_tau = inputs.multi_tau;  
    oracle = inputs.oracle;  


    %schedules
    prior_sched_prim = inputs.prior_sched_prim;
    prior_sched_sec = inputs.prior_sched_sec;
    trans_sched_prim = inputs.trans_sched_prim;
    trans_sched_sec = inputs.trans_sched_sec;

    %probes - prior
    exp_probes_features = inputs.exp_probes_features;
    exp_options_features = inputs.exp_options_features;
    obs_exp_choices = inputs.obs_exp_choices;
    
    inf_probes_features = inputs.inf_probes_features;
    inf_options_features = inputs.inf_options_features;
    obs_inf_choices = inputs.obs_inf_choices;
    
    %probes - transfer
    trans_exp_probes_features = inputs.trans_exp_probes_features;
    trans_exp_options_features = inputs.trans_exp_options_features;
    trans_obs_exp_choices = inputs.trans_obs_exp_choices;
    
    trans_inf_probes_features = inputs.trans_inf_probes_features;
    trans_inf_options_features = inputs.trans_inf_options_features;
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

    if transfer_model == 1
        w = params_to_opt(end); %mixing parameter
    else
        w = 0; %mixing parameter
    end

    all_subpch = [];
    k_all = zeros(dims_SF(1),1);

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
            sched_prim = prior_sched_prim;
            sched_sec = prior_sched_sec;

            all_exp_probes = exp_probes_features;
            all_exp_options = exp_options_features;
            all_exp_obs_choice = obs_exp_choices;

            all_inf_probes = inf_probes_features;
            all_inf_options = inf_options_features;
            all_inf_obs_choice = obs_inf_choices;

        else
            sched_prim = trans_sched_prim;
            sched_sec = trans_sched_sec;

            all_exp_probes = trans_exp_probes_features;
            all_exp_options = trans_exp_options_features;
            all_exp_obs_choice = trans_obs_exp_choices;

            all_inf_probes = trans_inf_probes_features;
            all_inf_options = trans_inf_options_features;
            all_inf_obs_choice = trans_obs_inf_choices;
        end

        %initialize successor representation matrix
        SF_M = zeros(dims_SF(1),dims_SF);
        sched = [reshape(sched_prim.',1,[])' reshape(sched_sec.',1,[])'];

        if multi_lr == 1            
            alpha_curr = alpha;
        else
            alpha_curr = alpha(task);
        end

        %% loop over trials/episodes
        iTrial = 0;
        k_all = zeros(size(SF_M,1),1);        
        
        for t = 1:(size(sched,1)-1)
            if t == 1
                SF_M(:,:,1) = SF_M + eye(dims_SF(1));
            end

            SF_M(:,:,t+1) = SF_M(:,:,t);

            twohot = zeros(1,dims_SF(1));
            twohot(sched(t,:)) = 1;

            for ii = 1:2
                
                FPE = .5 * alpha_curr * (twohot + gamma * sum(SF_M(sched(t+1,:),:,t),1) - sum(SF_M(sched(t,:),:,t),1));
                SF_M(sched(t,ii),:,t+1) = SF_M(sched(t,ii),:,t) + FPE;

                if task == 2
                    %count of state visitations, for exponential decay
                    k_all(sched(t,ii)) = k_all(sched(t,ii)) + 1;
                end
            end

            %% Probe trials
            if mod(t,size(sched_prim,2)) == 0 || t == size(sched,1)-1 %probes occur after every 12th trial, 
                % but we take the representation learned until then, not
                % considering the next transition (except last trial, which
                % happens at very end)

                if t == size(sched,1)-1 
                    SF_M_tmp = SF_M(:,:,end);
                else
                    SF_M_tmp = SF_M(:,:,end-1);
                end

                iTrial = iTrial + 1;

                if task == 2 && transfer_model == 1
                    if oracle == 0
                        %find optimal permutation of prior learning SR matrix
                        opt_perm = find_optimal_perm(M_prior, SF_M_tmp, 1000);
                        M_prior_permuted = M_prior(opt_perm,opt_perm);

                        % if iTrial == 36
                        % figure
                        % tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
                        % 
                        % nexttile
                        % imagesc(M_prior)
                        % colormap('hot')
                        % clim([0, 1.2]);
                        % yticks([]);
                        % yticklabels({});
                        % xticks([]);
                        % xticklabels({});
                        % pbaspect([2 2 1])
                        % 
                        % nexttile                        
                        % imagesc(M_prior_permuted)
                        % colormap('hot')
                        % clim([0, 1.2]);
                        % yticks([]);
                        % yticklabels({});
                        % xticks([]);
                        % xticklabels({});
                        % pbaspect([2 2 1])
                        % 
                        % nexttile                        
                        % imagesc(SF_M_tmp)
                        % colormap('hot')
                        % clim([0, 1.2]);
                        % yticks([]);
                        % yticklabels({});
                        % xticks([]);
                        % xticklabels({});
                        % pbaspect([2 2 1])
                        % set(gcf,'color','w');
                        % 
                        % exportgraphics(gcf, 'model_representations_features.png', 'Resolution', 1200);
                        % end

                    else
                        %Tryout mixing in actual TM of transfer learning phase to validate
                        if iTrial == 1
                            state_visit_counts = zeros(dims_SF);

                            if condition == 1

                                for i = 1:size(sched,1)-1

                                    current_state = sched(i,1);
                                    next_state = sched(i+1,1);

                                    %keep count of specific transition
                                    state_visit_counts(current_state,next_state) = state_visit_counts(current_state,next_state) + 1;
                                end

                            else

                                for i = 1:size(sched,1)-1

                                    current_state = sched(i,2);
                                    next_state = sched(i+1,2);

                                    %keep count of specific transition
                                    state_visit_counts(current_state,next_state) = state_visit_counts(current_state,next_state) + 1;
                                end

                            end

                            %make TM
                            TM = state_visit_counts ./ max(state_visit_counts,1);
                        end

                        M_prior_permuted = TM;
                    end
                    %% experience probes
                    %Softmax function arbitrating between alternative 1) next feature combination in the trajectory
                    %and alternative 2) some other feature combination in the trajectory
                    %calculate vector norm for presented options/successor features
                    succ_prob_prior(iTrial,:) = get_SF_probs(M_prior_permuted, all_exp_probes, all_exp_options, iTrial);
                    succ_prob_trans(iTrial,:) = get_SF_probs(SF_M_tmp, all_exp_probes, all_exp_options, iTrial);

                    %weighting of prior learning and transfer learning SF representation
                    %left option
                    succ_prob(iTrial,1) = w * succ_prob_prior(iTrial,1) + (1-w) * succ_prob_trans(iTrial,1);

                    %right option
                    succ_prob(iTrial,2) = w * succ_prob_prior(iTrial,2) + (1-w) * succ_prob_trans(iTrial,2);

                    %% inference probes
                    inf_succ_prob_prior(iTrial,:) = get_SF_probs(M_prior_permuted, all_inf_probes, all_inf_options, iTrial);
                    inf_succ_prob_trans(iTrial,:) = get_SF_probs(SF_M_tmp, all_inf_probes, all_inf_options, iTrial);

                    %weighting of prior learning and transfer learning SF representation
                    %left option
                    inf_succ_prob(iTrial,1) = w * inf_succ_prob_prior(iTrial,1) + (1-w) * inf_succ_prob_trans(iTrial,1);

                    %right option
                    inf_succ_prob(iTrial,2) = w * inf_succ_prob_prior(iTrial,2) + (1-w) * inf_succ_prob_trans(iTrial,2);

                else

                    succ_prob(iTrial,:) = get_SF_probs(SF_M_tmp, all_exp_probes, all_exp_options, iTrial);
                    inf_succ_prob(iTrial,:) = get_SF_probs(SF_M_tmp, all_inf_probes, all_inf_options, iTrial);

                end

                % if mod(iTrial,10) == 0
                %     if iTrial == 10
                %        figure;
                %        t = tiledlayout(1,4)
                %     end
                % 
                %     nexttile            
                %     imagesc(SF_M(:,:,end))
                % end

            end
        end
        


        %% Make choices
        for choice_type = 1:2 %experience or inference probe
            if multi_tau == 1            
                tau_curr = tau;
            else
                tau_curr = tau(choice_type);
            end

            % turn values into softmax choice probs
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
           
            if ~isempty(obs_choice)
                num_all_choices = size(all_subpch(~isnan(all_subpch)),1);
            else
                num_all_choices = [];
            end
        end

        %save prior learning SF matrix
        if task == 1
            M_prior = SF_M(:,:,end);
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

            if alpha < 0 || alpha > 1 || w < 0 || w > 1 || tau < 0 || tau > 100
                energy = 9999999;
            else
                energy = -nansum(log(all_subpch)); %calculate the negative log likelihood of parameters given the choices
            end

        else

            if alpha < 0 || alpha > 1 || w < 0 || w > 1 || tau(1) < 0 || tau(1) > 100 || tau(2) < 0 || tau(2) > 100
                energy = 9999999;
            else
                energy = -nansum(log(all_subpch)); %calculate the negative log likelihood of parameters given the choices
            end

        end

    elseif multi_lr == 2
        if multi_tau == 1
            
            if alpha(1) < 0 || alpha(1) > 1 || alpha(2) < 0 || alpha(2) > 1 || w < 0 || w > 1 || tau < 0 || tau > 100
                energy = 9999999;
            else
                energy = -nansum(log(all_subpch)); %calculate the negative log likelihood of parameters given the choices
            end

        else

            if alpha(1) < 0 || alpha(1) > 1 || alpha(2) < 0 || alpha(2) > 1 || w < 0 || w > 1 || tau(1) < 0 || tau(1) > 100 || tau(2) < 0 || tau(2) > 100
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

    %correct answers
    outputs.exp_probes_correct = inputs.exp_probes_correct;
    outputs.inf_probes_correct = inputs.inf_probes_correct;
    outputs.trans_exp_probes_correct = inputs.trans_exp_probes_correct;
    outputs.trans_inf_probes_correct = inputs.trans_inf_probes_correct;

    outputs.num_all_choices = num_all_choices;

