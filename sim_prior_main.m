%% Simulate different Successor Feature Models

clear all; close all; clc

format compact

addpath(genpath(['.', filesep, 'functions']))
addpath(genpath(['.', filesep, 'sim_parameters']))

%% settings
save_plots = 0;
save_results = 1;
n_subjects = 250; % number of subjects to be simulated
niter = 1; %repeat simulations how many times?

%fixing parameters
fixed_parameters = 0;
oracle = 0; %simulate "oracle model" (mixing in the true TM) in SF transfer model

% model_type = 5; %(1) SF compound model (compound-compund)
%(2) SF feature model (feature-feature)
%(5) SF feature model (feature-feature) with transfer learning , separate learning rate and mixing parameter

alpha_vec = [0.01 0.05]; %learning rate parameters used in simulation
tau_vec = [0.3 0.6]; %softmax stochasticity parameters used in simulation

plot_SR = 1; %if 1, will plot SF matrix plots across trials (increases run time), if 0, will not

%start parallel processing
delete(gcp('nocreate'))

num_workers = 8;
pool = parpool(num_workers);

for model_type = 1%[1 2 5]

    if model_type < 5
        transfer_model = 0;
    else
        transfer_model = 1;
    end

    for multi_lr = 2 %multiple (2 - prior/transfer) or a single learning rate parameter (1)
        for multi_tau = 1:2 %1 or 2 stochasticity parameters (overall or exp/inf probes separately)
            clear sim_res
            clear inf_sim_res
            clear sim_choices
            clear inf_sim_choices

            clear sim_res_trans
            clear inf_sim_res_trans
            clear sim_choices_trans
            clear inf_sim_choices_trans

            clear corr_resp_prim
            clear corr_resp_sec
            clear inf_corr_resp_prim
            clear inf_corr_resp_sec

            clear trans_corr_resp_prim
            clear trans_corr_resp_sec
            clear trans_inf_corr_resp_prim
            clear trans_inf_corr_resp_sec

            condition = randi([1 2],n_subjects,1);
            while sum(condition == 1) ~= sum(condition == 2)
                condition = randi([1 2],n_subjects,1);
            end

            res_dir = ['.', filesep, 'prior_simulations'];
            plotpath = ['.', filesep, 'prior_simulations', filesep, 'plots'];

            %% parameters
            gamma = .9; %discounting factor gamma in SF update, less updating of more distal states (less likely to visit states further away)

            %load subject-specific parameters
            space_alpha_ini = load('space_alpha.mat');
            space_alpha_full = space_alpha_ini.space_alpha;
            space_alpha_full_unsorted = space_alpha_ini.space_alpha;
            [~, param_idx] = sort(space_alpha_full);
            space_alpha_full = space_alpha_full(param_idx,:);

            if fixed_parameters == 1
                space_alpha_full = repmat(0.12,size(space_alpha_full_unsorted));
            end


            space_tau_ini = load('space_tau.mat');
            space_tau_full = space_tau_ini.space_tau;
            space_tau_full_unsorted = space_tau_ini.space_tau;
            space_tau_full = space_tau_full(param_idx,:);


            space_w_ini = load('space_omega.mat');
            space_w_full = space_w_ini.space_omega;
            space_w_full_unsorted = space_w_ini.space_omega;
            space_w_full = space_w_full(param_idx,:);

            if fixed_parameters == 2
                space_w_full = repmat(0.12,size(space_w_full_unsorted));
            end


            space_omega_ini = load('space_omega.mat');
            space_omega_full = space_omega_ini.space_omega;
            space_omega_full_unsorted = space_omega_ini.space_omega;
            space_omega_full = space_omega_full(param_idx,:);

            %select randomly drawn number of parameters (number of subjects in the
            %simulation)
            rand_idx = randperm(size(space_alpha_full,1),n_subjects);
            space_alpha = zeros(size(space_alpha_full,1)/2,1);
            space_tau = zeros(size(space_tau_full,1)/2,1);
            space_w = zeros(size(space_w_full,1)/2,1);
            space_omega = zeros(size(space_omega_full,1)/2,1);

            %match parameters across conditions
            space_alpha_half = space_alpha_full_unsorted(1:size(space_alpha_full,1)/2);
            space_tau_half = space_tau_full_unsorted(1:size(space_alpha_full,1)/2);
            space_w_half = space_w_full_unsorted(1:size(space_w_full,1)/2);
            space_omega_half = space_omega_full_unsorted(1:size(space_omega_full,1)/2);

            %condition 1
            for k = 1:size(condition,1)
                if condition(k) == 1
                    space_alpha(k) = space_alpha_half(1);
                    space_tau(k) = space_tau_half(1);
                    space_w(k) = space_w_half(1);
                    space_omega(k) = space_omega_half(1);

                    space_alpha_half(1) = [];
                    space_tau_half(1) = [];
                    space_w_half(1) = [];
                    space_omega_half(1) = [];
                end

            end

            space_alpha_half = space_alpha_full_unsorted(1:size(space_alpha_full,1)/2);
            space_tau_half = space_tau_full_unsorted(1:size(space_alpha_full,1)/2);
            space_w_half = space_w_full_unsorted(1:size(space_w_full,1)/2);
            space_omega_half = space_omega_full_unsorted(1:size(space_omega_full,1)/2);

            %condition 2
            for k = 1:size(condition,1)
                if condition(k) == 2
                    space_alpha(k) = space_alpha_half(1);
                    space_tau(k) = space_tau_half(1);
                    space_w(k) = space_w_half(1);
                    space_omega(k) = space_omega_half(1);

                    space_alpha_half(1) = [];
                    space_tau_half(1) = [];
                    space_w_half(1) = [];
                    space_omega_half(1) = [];
                end

            end

            space_w = repmat(0.25, size(space_w,1),1);

            if multi_lr == 2
                space_alpha = repmat(alpha_vec, size(space_alpha,1),1);
            else
                space_alpha = repmat(alpha_vec(1), size(space_alpha,1),1);
            end

            if multi_tau == 2
                space_tau = repmat(tau_vec, size(space_tau,1),1);
            else
                space_tau = repmat(tau_vec(1), size(space_tau,1),1);
            end

            %nodes of the factorized graph in learning phase
            feature_names = {'A', 'B', 'C', 'D', ...
                '1', '2', '3', '4', '5', '6'};

            %load trials/schedule
            cd(['.', filesep, 'schedules'])

            all_scheds = dir('Full*');

            cd ..
            %load schedule number
            load('sched_num_all.mat')

            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Model simulation

            for subj = 1:length(space_alpha)
                alpha = space_alpha(subj,:); %learning rate(s)
                tau = space_tau(subj,:); %Softmax stochasticity parameter
                w_mix = space_w(subj); %mixing parameter
                omega = space_omega(subj);

                % subj = 7; %schedule where transfer is easy to see (intuitive
                % permutation of states)
                % condition(subj) = 1;

                % load schedules
                if condition(subj) == 1
                    sched_num = sched_num_all(subj);
                else
                    sched_num = sched_num_all(subj) + 8;
                end

                cd(['.', filesep, 'schedules'])
                sched_raw = readtable(all_scheds(sched_num).name);
                cd ..

                task_begin_idx = repmat(1,length(space_alpha),1);
                task_end_idx  = repmat(size(sched_raw,1),length(space_alpha),1);

                sched_raw.ParticipantPrivateID = repmat(6577523,size(sched_raw,1),1);

                %% Preprocessing/preparation of data
                preproc_outputs = preprocess_data(sched_raw, task_begin_idx, task_end_idx, subj, feature_names, model_type);

                if model_type > 1
                    %check if choice trials are correct
                    correct_trials = check_choice_trials(preproc_outputs, model_type);

                    all_choice_trials_correct(subj,1) = sum(correct_trials.prior_exp);
                    all_choice_trials_correct(subj,2) = sum(correct_trials.prior_inf);
                    all_choice_trials_correct(subj,3) = sum(correct_trials.trans_exp);
                    all_choice_trials_correct(subj,4) = sum(correct_trials.trans_inf);

                    if subj == length(space_alpha)
                        sum(all_choice_trials_correct(:,1) == 48) == n_subjects
                        sum(all_choice_trials_correct(:,2) == 48) == n_subjects
                        sum(all_choice_trials_correct(:,3) == 36) == n_subjects
                        sum(all_choice_trials_correct(:,4) == 36) == n_subjects
                    end
                end

                %prepare model inputs
                preproc_outputs.simulation = 1;
                preproc_outputs.fixed_params = fixed_parameters;
                preproc_outputs.gamma = gamma;
                preproc_outputs.multi_lr = multi_lr;
                preproc_outputs.multi_tau = multi_tau;
                preproc_outputs.condition = condition(subj);
                preproc_outputs.oracle = oracle;
                preproc_outputs.transfer_model = transfer_model;


                % parfor iter = 1:niter
                    for iter = 1:niter


                    %% SF compound model
                    if model_type == 1

                        params = [alpha tau];

                        %run simulation
                        [~, model_outputs] = SF_compound_model_all(preproc_outputs, params);


                        %% SF feature model
                    elseif model_type == 2

                        params = [alpha tau];


                        %run simulation
                        [~, model_outputs] = SF_feature_model_all(preproc_outputs, params);

                        %% SF feature model with transfer learning, separate learning rate and mixing parameter
                    elseif model_type == 5

                        params = [alpha tau w_mix];

                        % other_params = [0.1795    0.6575    0.9000];
                        % params = [0.1    0.3  0.001];

                        %run simulation
                        [~, model_outputs] = SF_feature_model_all(preproc_outputs, params);

                    end

                    %% store trial type specific simulated choice and prob correct

                    sim_res(:,subj,iter) = model_outputs.prior_sim_res_corr_exp;
                    inf_sim_res(:,subj,iter) = model_outputs.prior_sim_res_corr_inf;
                    sim_choices(:,subj,iter) = model_outputs.prior_sim_choices_exp;
                    inf_sim_choices(:,subj,iter) =  model_outputs.prior_sim_choices_inf;

                    sim_res_trans(:,subj,iter) = model_outputs.trans_sim_res_corr_exp;
                    inf_sim_res_trans(:,subj,iter) = model_outputs.trans_sim_res_corr_inf;
                    sim_choices_trans(:,subj,iter) = model_outputs.trans_sim_choices_exp;
                    inf_sim_choices_trans(:,subj,iter) = model_outputs.trans_sim_choices_inf;


                    corr_resp_prim(subj,:,iter) = model_outputs.prior_sim_res_corr_exp(preproc_outputs.exp_trial_type==1);
                    corr_resp_sec(subj,:,iter) = model_outputs.prior_sim_res_corr_exp(preproc_outputs.exp_trial_type==2);

                    inf_corr_resp_prim(subj,:,iter) = model_outputs.prior_sim_res_corr_inf(preproc_outputs.inf_trial_type==1);
                    inf_corr_resp_sec(subj,:,iter) = model_outputs.prior_sim_res_corr_inf(preproc_outputs.inf_trial_type==2);


                    trans_corr_resp_prim(subj,:,iter) = model_outputs.trans_sim_res_corr_exp(preproc_outputs.trans_exp_trial_type==1);
                    trans_corr_resp_sec(subj,:,iter) = model_outputs.trans_sim_res_corr_exp(preproc_outputs.trans_exp_trial_type==2);

                    trans_inf_corr_resp_prim(subj,:,iter) = model_outputs.trans_sim_res_corr_inf(preproc_outputs.trans_inf_trial_type==1);
                    trans_inf_corr_resp_sec(subj,:,iter) = model_outputs.trans_sim_res_corr_inf(preproc_outputs.trans_inf_trial_type==2);



                end %end iteration parfor loop


                % if model_type > 1
                % %% plot predicted probability
                % for j = 1:size(model_outputs.prior_chprob_exp,1)
                %     pp_exp_prior(j,1) = model_outputs.prior_chprob_exp(j,model_outputs.exp_probes_correct(j,1));
                %     pp_inf_prior(j,1) = model_outputs.prior_chprob_inf(j,model_outputs.inf_probes_correct(j,1));
                % end
                %
                % for j = 1:size(model_outputs.trans_chprob_exp,1)
                %     pp_exp_trans(j,1) = model_outputs.trans_chprob_exp(j,model_outputs.trans_exp_probes_correct(j,1));
                %     pp_inf_trans(j,1) = model_outputs.trans_chprob_inf(j,model_outputs.trans_inf_probes_correct(j,1));
                % end


                % figure
                % subplot(2,2,1)
                % plot(pp_exp_prior(preproc_outputs.exp_trial_type == 1), 'g'); hold on
                % plot(pp_exp_prior(preproc_outputs.exp_trial_type == 2), 'k')
                %
                % subplot(2,2,2)
                % plot(pp_inf_prior(preproc_outputs.inf_trial_type == 1), 'g'); hold on
                % plot(pp_inf_prior(preproc_outputs.inf_trial_type == 2), 'k')
                %
                % subplot(2,2,3)
                % plot(pp_exp_trans(preproc_outputs.trans_exp_trial_type == 1), 'g'); hold on
                % plot(pp_exp_trans(preproc_outputs.trans_exp_trial_type == 2), 'k')
                %
                % subplot(2,2,4)
                % plot(pp_inf_trans(preproc_outputs.trans_inf_trial_type == 1), 'g'); hold on
                % plot(pp_inf_trans(preproc_outputs.trans_inf_trial_type == 2), 'k')
                %
                % % suptitle(num2str(params))
                % sgtitle(['condition: ' num2str(condition(subj)) ', params: ' num2str(params)]);
                % end





            end %end subject loop

            mean(space_alpha(condition == 1))
            std(space_alpha(condition == 1))

            mean(space_alpha(condition == 2))
            std(space_alpha(condition == 2))

            mean(space_w(condition == 1))
            std(space_w(condition == 1))

            mean(space_w(condition == 2))
            std(space_w(condition == 2))



            %compute averages
            sim_res = mean(sim_res,3);
            inf_sim_res = mean(inf_sim_res,3);
            sim_res_trans = mean(sim_res_trans,3);
            inf_sim_res_trans = mean(inf_sim_res_trans,3);

            corr_resp_prim = mean(corr_resp_prim,3);
            corr_resp_sec = mean(corr_resp_sec,3);
            inf_corr_resp_prim = mean(inf_corr_resp_prim,3);
            inf_corr_resp_sec = mean(inf_corr_resp_sec,3);

            trans_corr_resp_prim = mean(trans_corr_resp_prim,3);
            trans_corr_resp_sec = mean(trans_corr_resp_sec,3);

            trans_inf_corr_resp_prim = mean(trans_inf_corr_resp_prim,3);
            trans_inf_corr_resp_sec = mean(trans_inf_corr_resp_sec,3);


            %% Plots

            %% plot percent correct choices across probe trials - learning trajectories
            %some example subjects
            avg_acc = mean(sim_res,1);
            [~,ranking]=ismember(avg_acc,sort(avg_acc,'descend'));
            blub = [avg_acc; ranking];
            sim_res_sorted = [sim_res' ranking'];
            sim_res_sorted = sortrows(sim_res_sorted,size(sim_res_sorted,2));
            sim_res_sorted = sim_res_sorted';

            traject_to_plot = sim_res_sorted(:, [1 round(n_subjects/3.75) n_subjects/2 n_subjects-(round(n_subjects/3.75)) n_subjects]);

            names_alpha = space_alpha([1 round(n_subjects/3.75) n_subjects/2 n_subjects-(round(n_subjects/3.75)) n_subjects]);
            names_tau = space_tau([1 round(n_subjects/3.75) n_subjects/2 n_subjects-(round(n_subjects/3.75)) n_subjects]);

            mean(avg_acc(condition == 1))
            mean(avg_acc(condition == 2))

            %average plot across simulations
            corr_resp_all = [NaN(size(sim_res,2),1) sim_res'];
            inf_corr_resp_all = [NaN(size(inf_sim_res,2),1) inf_sim_res'];

            cols = [];

            cols.k = [0 0 0];
            cols.b = [0 .058 .686];
            cols.y = [1  .828 0];
            cols.grey = [0.7843 0.7843 0.7843];
            cols.dgrey = [0.1922 0.2000 0.2078];
            cols.green = [0 186 85]/255;
            cols.lblue = [10 118 194]/255;
            cols.ytrend = [150 150 105]/255;
            cols.btrend = [110 110 150]/255;
            cols.gtrend = [80 120 90]/255;


            % Fit linear trends to data
            %overall SL - prior
            linearCoefficients_sl = polyfit(1:48,nanmean(corr_resp_all(:,2:end)),1); % First degree linear fit
            yFit_sl_prior = polyval(linearCoefficients_sl, 1:48);
            %overall inf - prior
            linearCoefficients_inf = polyfit(1:48,nanmean(inf_corr_resp_all(:,2:end)),1); % First degree linear fit
            yFit_inf_prior = polyval(linearCoefficients_inf, 1:48);

            h = figure('visible','on');
            %subplot(1,3,2)
            subplot(1,2,1)
            stdshade(corr_resp_all(:,2:end), .4, cols.lblue, [], 5); hold on;
            stdshade(inf_corr_resp_all(:,2:end), .4, cols.green, [], 5);
            set(findobj(gca,'type','line'),'linew',2);
            plot(1:48,yFit_sl_prior, '--','Color', cols.btrend)
            plot(1:48,yFit_inf_prior, '--','Color', cols.gtrend)
            hold off
            box off
            xbounds = [0 size(corr_resp_all,2)-1];
            ylim([0.25 1]);
            set(findobj(gca,'type','line'),'linew',3)
            set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
            set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
            set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
            yline(0.5,'k--');
            set(gcf,'color','w');
            set(gca,'TickDir','out')
            ylabel('Probability of Correct Answers', 'FontSize',20,...
                'Color','k')
            title('Average correct responses')
            xtickangle(45)
            legend('','Exp Probes', '','Inf Probes','Location','NorthWestOutside', 'FontSize',12)



            %%average p_corr for trial types
            corr_resp_prim = [NaN(size(corr_resp_prim,1),1) corr_resp_prim];
            corr_resp_sec = [NaN(size(corr_resp_sec,1),1) corr_resp_sec];

            % Fit linear trends to data - SL
            %overall primary
            linearCoefficients_prim_sl = polyfit(1:24,nanmean(corr_resp_prim(:,2:end)),1); % First degree linear fit
            yFit_prim_sl = polyval(linearCoefficients_prim_sl, 1:24);

            %overall secondary
            linearCoefficients_sec_sl = polyfit(1:24,nanmean(corr_resp_sec(:,2:end)),1); % First degree linear fit
            yFit_sec_sl = polyval(linearCoefficients_sec_sl, 1:24);

            subplot(2,2,2)
            stdshade(corr_resp_prim(:,2:end), .1, cols.y, [], 1); hold on;
            stdshade(corr_resp_sec(:,2:end), .1, cols.b, [], 1);
            set(findobj(gca,'type','line'),'linew',2);
            plot(1:24,yFit_prim_sl, '--','Color', cols.ytrend)
            plot(1:24,yFit_sec_sl, '--','Color', cols.btrend)
            hold off
            box off
            xbounds = [0 size(corr_resp_prim,2)-1];
            ylim([0.25 1]);
            set(findobj(gca,'type','line'),'linew',3)
            set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
            set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
            set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
            yline(0.5,'k--');
            set(gcf,'color','w');
            set(gca,'TickDir','out')
            % ylabel('Probability of Correct Answers', 'FontSize',20,...
            %     'Color','k')
            title('Correct responses by trial type')
            xtickangle(45)
            title('Exp Probes: Trial type')
            legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)


            inf_corr_resp_prim = [NaN(size(inf_corr_resp_prim,1),1) inf_corr_resp_prim];
            inf_corr_resp_sec = [NaN(size(inf_corr_resp_sec,1),1) inf_corr_resp_sec];

            % Fit linear trends to data - inference
            %overall primary
            linearCoefficients_prim = polyfit(1:24,nanmean(inf_corr_resp_prim(:,2:end)),1); % First degree linear fit
            yFit_prim = polyval(linearCoefficients_prim, 1:24);

            %overall secondary
            linearCoefficients_sec = polyfit(1:24,nanmean(inf_corr_resp_sec(:,2:end)),1); % First degree linear fit
            yFit_sec = polyval(linearCoefficients_sec, 1:24);

            %%average p_corr for trial types - inference
            subplot(2,2,4)
            stdshade(inf_corr_resp_prim(:,2:size(inf_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
            stdshade(inf_corr_resp_sec(:,2:size(inf_corr_resp_sec,2)), .1, cols.b, [], 1);
            set(findobj(gca,'type','line'),'linew',2);
            plot(1:24,yFit_prim, '--','Color', cols.ytrend)
            plot(1:24,yFit_sec, '--','Color', cols.btrend)
            hold off
            box off
            xbounds = [0 size(inf_corr_resp_prim,2)-1];
            ylim([0.25 1]);
            set(findobj(gca,'type','line'),'linew',3)
            set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
            set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
            set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
            yline(0.5,'k--');
            set(gcf,'color','w');
            set(gca,'TickDir','out')
            % ylabel('Probability of Correct Answers', 'FontSize',20,...
            %     'Color','k')
            title('Inference Probes: Trial type')
            legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
            xtickangle(45)

            if save_plots == 1
                screen_size = get(0, 'ScreenSize');
                origSize = get(h, 'Position'); % grab original on screen size
                set(h, 'Position', [0 0 screen_size(3) screen_size(4)/2 ] ); %set to scren size
                set(h,'PaperPositionMode','auto') %set paper pos for printing
                % saveas(h, 'fig_trajectories.png') % save figure
                ax = h;
                exportgraphics(ax,'fig_trajectories.png','Resolution',1200)
                set(h,'Position', origSize) %set back to original dimensions
            end



            %% Transfer
            % Trajectories in transfer overall - split up for conditions
            trial_num_transfer = 36;

            %average plot across simulations
            trans_corr_resp_all = [NaN(size(sim_res_trans,2),1) sim_res_trans'];
            trans_inf_corr_resp_all = [NaN(size(inf_sim_res_trans,2),1) inf_sim_res_trans'];

            trans_corr_resp_prim = [NaN(size(trans_corr_resp_prim,1),1) trans_corr_resp_prim];
            trans_corr_resp_sec = [NaN(size(trans_corr_resp_sec,1),1) trans_corr_resp_sec];

            trans_inf_corr_resp_prim = [NaN(size(trans_inf_corr_resp_prim,1),1) trans_inf_corr_resp_prim];
            trans_inf_corr_resp_sec = [NaN(size(trans_inf_corr_resp_sec,1),1) trans_inf_corr_resp_sec];

            %average p_corr - SL and inference probes - transfer learning

            % Fit linear trends to data
            %overall SL
            linearCoefficients_trans_sl_cond1 = polyfit(1:trial_num_transfer,nanmean(trans_corr_resp_all(condition == 1,2:end)),1); % First degree linear fit
            yFit_trans_sl_cond1 = polyval(linearCoefficients_trans_sl_cond1, 1:trial_num_transfer);

            linearCoefficients_trans_sl_cond2 = polyfit(1:trial_num_transfer,nanmean(trans_corr_resp_all(condition == 2,2:end)),1); % First degree linear fit
            yFit_trans_sl_cond2 = polyval(linearCoefficients_trans_sl_cond2, 1:trial_num_transfer);

            %overall inf
            linearCoefficients_trans_inf_cond1 = polyfit(1:trial_num_transfer,nanmean(trans_inf_corr_resp_all(condition == 1,2:end)),1); % First degree linear fit
            yFit_trans_inf_cond1 = polyval(linearCoefficients_trans_inf_cond1, 1:trial_num_transfer);

            linearCoefficients_trans_inf_cond2 = polyfit(1:trial_num_transfer,nanmean(trans_inf_corr_resp_all(condition == 2,2:end)),1); % First degree linear fit
            yFit_trans_inf_cond2 = polyval(linearCoefficients_trans_inf_cond2, 1:trial_num_transfer);


            h = figure('visible','on');
            %subplot(1,3,2)
            subplot(1,2,1)
            stdshade(trans_corr_resp_all(condition == 1,2:size(trans_corr_resp_all,2)), .4, cols.lblue, [], 5); hold on;
            stdshade(trans_inf_corr_resp_all(condition == 1,2:size(trans_inf_corr_resp_all,2)), .4, cols.green, [], 5);
            set(findobj(gca,'type','line'),'linew',2);
            plot(1:trial_num_transfer,yFit_trans_sl_cond1, '--','Color', cols.btrend)
            plot(1:trial_num_transfer,yFit_trans_inf_cond1, '--','Color', cols.gtrend)
            hold off
            box off
            xbounds = [0 size(trans_corr_resp_all,2)-1];
            ylim([0.25 1]);
            xlim(xbounds);
            set(findobj(gca,'type','line'),'linew',3)
            set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
            set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
            set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
            yline(0.5,'k--');
            set(gcf,'color','w');
            set(gca,'TickDir','out')
            ylabel('Probability of Correct Answers', 'FontSize',20,...
                'Color','k')
            title('Average correct responses - Transfer, cond1')
            xtickangle(45)
            %legend('','Seq. Probes', '','Inf Probes','Location','NorthWestOutside', 'FontSize',12)

            subplot(1,2,2)
            stdshade(trans_corr_resp_all(condition == 2,2:size(trans_corr_resp_all,2)), .4, cols.lblue, [], 5); hold on;
            stdshade(trans_inf_corr_resp_all(condition == 2,2:size(trans_inf_corr_resp_all,2)), .4, cols.green, [], 5);
            set(findobj(gca,'type','line'),'linew',2);
            plot(1:trial_num_transfer,yFit_trans_sl_cond2, '--','Color', cols.btrend)
            plot(1:trial_num_transfer,yFit_trans_inf_cond2, '--','Color', cols.gtrend)
            hold off
            box off
            xbounds = [0 size(trans_corr_resp_all,2)-1];
            ylim([0.25 1]);
            xlim(xbounds);
            set(findobj(gca,'type','line'),'linew',3)
            set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
            set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
            set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
            yline(0.5,'k--');
            set(gcf,'color','w');
            set(gca,'TickDir','out')
            ylabel('Probability of Correct Answers', 'FontSize',20,...
                'Color','k')
            title('Average correct responses - Transfer, cond2')
            xtickangle(45)

            if save_plots == 1
                screen_size = get(0, 'ScreenSize');
                origSize = get(h, 'Position'); % grab original on screen size
                set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size
                set(h,'PaperPositionMode','auto') %set paper pos for printing
                % saveas(h, 'fig_trajectories.png') % save figure
                ax = h;
                exportgraphics(ax,'fig_trajectories_transfer.png','Resolution',1200)
                set(h,'Position', origSize) %set back to original dimensions
            end



            %% Trajectories in transfer - primary/secondary - split up for conditions
            % Fit linear trends to data - SL
            %overall primary
            linearCoefficients_prim_sl_trans_cond1 = polyfit(1:trial_num_transfer/2,nanmean(trans_corr_resp_prim(condition == 1,2:end)),1); % First degree linear fit
            yFit_prim_sl_trans_cond1 = polyval(linearCoefficients_prim_sl_trans_cond1, 1:trial_num_transfer/2);

            [linearCoefficients_prim_sl_trans_cond2, S] = polyfit(1:trial_num_transfer/2,nanmean(trans_corr_resp_prim(condition == 2,2:end)),1); % First degree linear fit
            yFit_prim_sl_trans_cond2 = polyval(linearCoefficients_prim_sl_trans_cond2, 1:trial_num_transfer/2);
            CI = polyparci(linearCoefficients_prim_sl_trans_cond2,S);

            %overall secondary
            linearCoefficients_sec_sl_trans_cond1 = polyfit(1:trial_num_transfer/2,nanmean(trans_corr_resp_sec(condition == 1,2:end)),1); % First degree linear fit
            yFit_sec_sl_trans_cond1 = polyval(linearCoefficients_sec_sl_trans_cond1, 1:trial_num_transfer/2);

            linearCoefficients_sec_sl_trans_cond2 = polyfit(1:trial_num_transfer/2,nanmean(trans_corr_resp_sec(condition == 2,2:end)),1); % First degree linear fit
            yFit_sec_sl_trans_cond2 = polyval(linearCoefficients_sec_sl_trans_cond2, 1:trial_num_transfer/2);

            %%average p_corr for trial types
            %condition 1
            h = figure('visible','on');
            subplot(2,2,1)
            stdshade(trans_corr_resp_prim(condition == 1,2:size(trans_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
            stdshade(trans_corr_resp_sec(condition == 1,2:size(trans_corr_resp_sec,2)), .1, cols.b, [], 1);
            set(findobj(gca,'type','line'),'linew',2);
            plot(1:trial_num_transfer/2,yFit_prim_sl_trans_cond1, '--','Color', cols.ytrend)
            plot(1:trial_num_transfer/2,yFit_sec_sl_trans_cond1, '--','Color', cols.btrend)
            hold off
            box off
            xbounds = [0 size(trans_corr_resp_prim,2)-1];
            ylim([0.25 1]);
            xlim(xbounds);
            set(findobj(gca,'type','line'),'linew',3)
            set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
            set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
            set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
            yline(0.5,'k--');
            set(gcf,'color','w');
            set(gca,'TickDir','out')
            % ylabel('Probability of Correct Answers', 'FontSize',20,...
            %     'Color','k')
            title('Exp probes: Trial type - Transfer, cond1')
            legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
            xtickangle(45)


            %%average p_corr for trial types
            %condition 2
            subplot(2,2,3)
            stdshade(trans_corr_resp_prim(condition == 2,2:size(trans_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
            stdshade(trans_corr_resp_sec(condition == 2,2:size(trans_corr_resp_sec,2)), .1, cols.b, [], 1);
            set(findobj(gca,'type','line'),'linew',2);
            plot(1:trial_num_transfer/2,yFit_prim_sl_trans_cond2, '--','Color', cols.ytrend)
            plot(1:trial_num_transfer/2,yFit_sec_sl_trans_cond2, '--','Color', cols.btrend)
            hold off
            box off
            xbounds = [0 size(trans_corr_resp_prim,2)-1];
            ylim([0.25 1]);
            xlim(xbounds);
            set(findobj(gca,'type','line'),'linew',3)
            set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
            set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
            set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
            yline(0.5,'k--');
            set(gcf,'color','w');
            set(gca,'TickDir','out')
            % ylabel('Probability of Correct Answers', 'FontSize',20,...
            %     'Color','k')
            title('Exp probes: Trial type - Transfer, cond2')
            legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
            xtickangle(45)



            %%average p_corr for trial types - inference, transfer
            %overall primary
            linearCoefficients_prim_inf_trans_cond1 = polyfit(1:trial_num_transfer/2,nanmean(trans_inf_corr_resp_prim(condition == 1,2:end)),1); % First degree linear fit
            yFit_prim_inf_trans_cond1 = polyval(linearCoefficients_prim_inf_trans_cond1, 1:trial_num_transfer/2);

            linearCoefficients_prim_inf_trans_cond2 = polyfit(1:trial_num_transfer/2,nanmean(trans_inf_corr_resp_prim(condition == 2,2:end)),1); % First degree linear fit
            yFit_prim_inf_trans_cond2 = polyval(linearCoefficients_prim_inf_trans_cond2, 1:trial_num_transfer/2);

            %overall secondary
            linearCoefficients_sec_inf_trans_cond1 = polyfit(1:trial_num_transfer/2,nanmean(trans_inf_corr_resp_sec(condition == 1,2:end)),1); % First degree linear fit
            yFit_sec_inf_trans_cond1 = polyval(linearCoefficients_sec_inf_trans_cond1, 1:trial_num_transfer/2);

            linearCoefficients_sec_inf_trans_cond2 = polyfit(1:trial_num_transfer/2,nanmean(trans_inf_corr_resp_sec(condition == 2,2:end)),1); % First degree linear fit
            yFit_sec_inf_trans_cond2 = polyval(linearCoefficients_sec_inf_trans_cond2, 1:trial_num_transfer/2);

            %condition 1
            subplot(2,2,2)
            stdshade(trans_inf_corr_resp_prim(condition == 1,2:size(trans_inf_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
            stdshade(trans_inf_corr_resp_sec(condition == 1,2:size(trans_inf_corr_resp_sec,2)), .1, cols.b, [], 1);
            set(findobj(gca,'type','line'),'linew',2);
            plot(1:trial_num_transfer/2,yFit_prim_inf_trans_cond1, '--','Color', cols.ytrend)
            plot(1:trial_num_transfer/2,yFit_sec_inf_trans_cond1, '--','Color', cols.btrend)
            hold off
            box off
            xbounds = [0 size(trans_inf_corr_resp_prim,2)-1];
            ylim([0.25 1]);
            xlim(xbounds);
            set(findobj(gca,'type','line'),'linew',3)
            set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
            set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
            set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
            yline(0.5,'k--');
            set(gcf,'color','w');
            set(gca,'TickDir','out')
            % ylabel('Probability of Correct Answers', 'FontSize',20,...
            %     'Color','k')
            title('Inference probes: Trial type - Transfer, cond1')
            legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
            xtickangle(45)

            %condition 2
            subplot(2,2,4)
            stdshade(trans_inf_corr_resp_prim(condition == 2,2:size(trans_inf_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
            stdshade(trans_inf_corr_resp_sec(condition == 2,2:size(trans_inf_corr_resp_sec,2)), .1, cols.b, [], 1);
            set(findobj(gca,'type','line'),'linew',2);
            plot(1:trial_num_transfer/2,yFit_prim_inf_trans_cond2, '--','Color', cols.ytrend)
            plot(1:trial_num_transfer/2,yFit_sec_inf_trans_cond2, '--','Color', cols.btrend)
            hold off
            box off
            xbounds = [0 size(trans_inf_corr_resp_prim,2)-1];
            ylim([0.25 1]);
            xlim(xbounds);
            set(findobj(gca,'type','line'),'linew',3)
            set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
            set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
            set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
            yline(0.5,'k--');
            set(gcf,'color','w');
            set(gca,'TickDir','out')
            % ylabel('Probability of Correct Answers', 'FontSize',20,...
            %     'Color','k')
            title('Inference probes: Trial type - Transfer, cond2')
            legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
            xtickangle(45)

            if save_plots == 1
                screen_size = get(0, 'ScreenSize');
                origSize = get(h, 'Position'); % grab original on screen size
                set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size
                set(h,'PaperPositionMode','auto') %set paper pos for printing
                % saveas(h, 'fig_trajectories.png') % save figure
                ax = h;
                exportgraphics(ax,'fig_trajectories_trial_type_transfer.png','Resolution',1200)
                set(h,'Position', origSize) %set back to original dimensions
            end

            %% bar plots
            barColorMap = [cols.grey; cols.dgrey; cols.b; cols.y];

            %averages in different trial types
            %experience probes
            %prior
            avg_acc_primary = mean(corr_resp_prim(:,2:end),2);
            avg_acc_secondary = mean(corr_resp_sec(:,2:end),2);

            %transfer
            trans_avg_acc = mean(trans_corr_resp_all(:,2:end),2);
            trans_avg_acc_primary = mean(trans_corr_resp_prim(:,2:end),2);
            trans_avg_acc_secondary = mean(trans_corr_resp_sec(:,2:end),2);

            exp_res_to_plot = [avg_acc' avg_acc_primary avg_acc_secondary trans_avg_acc trans_avg_acc_primary trans_avg_acc_secondary];

            % inference probes
            %prior
            avg_acc_inf_primary = mean(inf_corr_resp_prim(:,2:end),2);
            avg_acc_inf_secondary = mean(inf_corr_resp_sec(:,2:end),2);
            avg_acc_inf = mean(inf_sim_res,1);

            %transfer
            trans_avg_acc_inf_primary = mean(trans_inf_corr_resp_prim(:,2:end),2);
            trans_avg_acc_inf_secondary = mean(trans_inf_corr_resp_sec(:,2:end),2);
            trans_avg_acc_inf = mean(inf_sim_res_trans,1);

            res_to_plot = [avg_acc' avg_acc_inf' trans_avg_acc trans_avg_acc_inf'];
            inf_to_plot = [avg_acc_inf' avg_acc_inf_primary avg_acc_inf_secondary trans_avg_acc_inf' trans_avg_acc_inf_primary trans_avg_acc_inf_secondary];


            %% Plotting
            ylim_lower = 0.4;
            y_ticks = ylim_lower:.2:1;

            results(:,8) = condition;

            plotting_simulations

            %save parameter estimates
            if multi_lr == 1
                if multi_tau == 1
                    mod_name = '1Alpha_1Tau';
                else
                    mod_name = '1Alpha_2Tau';
                end
            elseif multi_lr == 2
                if multi_tau == 1
                    mod_name = '2Alpha_1Tau';
                else
                    mod_name = '2Alpha_2Tau';
                end
            end

            if model_type == 1
                model_name = ['SF_compounds_', num2str(mod_name), '_sim_data.mat'];

            elseif model_type == 2
                model_name = ['SF_features_', num2str(mod_name), '_sim_data.mat'];

            elseif model_type == 5
                model_name = ['SF_features_transfer_', num2str(mod_name), '_sim_data.mat'];

            end

            if save_plots == 1
                screen_size = get(0, 'ScreenSize');
                origSize = get(h, 'Position'); % grab original on screen size
                set(h, 'Position', [0 0 screen_size(3) screen_size(4)/2 ] ); %set to scren size
                set(h,'PaperPositionMode','auto') %set paper pos for printing
                % saveas(h, 'fig_trajectories.png') % save figure
                ax = h;
                exportgraphics(ax,[plotpath, filesep, 'fig_pcorrect_all_' model_name '.png'],'Resolution',300)
                set(h,'Position', origSize) %set back to original dimensions
            end

            %% heatmap of performance, per condition
            figure

            subplot(2,2,1)
            acc_heatmap_plot(exp_res_to_plot, results, 1, 1);

            subplot(2,2,2)
            acc_heatmap_plot(inf_to_plot, results, 1, 2);

            subplot(2,2,3)
            acc_heatmap_plot(exp_res_to_plot, results, 2, 1);

            subplot(2,2,4)
            acc_heatmap_plot(inf_to_plot, results, 2, 2);

            cd(plotpath)
            if save_plots == 1
                screen_size = get(0, 'ScreenSize');
                origSize = get(h, 'Position'); % grab original on screen size
                set(h, 'Position', [0 0 screen_size(3)/2 screen_size(4)/1.5 ] ); %set to scren size
                set(h,'PaperPositionMode','auto') %set paper pos for printing
                % saveas(h, 'fig_trajectories.png') % save figure
                ax = h;
                exportgraphics(ax,['heatmap_prior_simulation_model_' num2str(model_type) '.png'],'Resolution',600)
                set(h,'Position', origSize) %set back to original dimensions
            end

            cd ..

            %% save data for plotting
            sim_data.results = results;
            sim_data.res_to_plot = res_to_plot;
            sim_data.exp_res_to_plot = exp_res_to_plot;
            sim_data.inf_to_plot = inf_to_plot;

            %save parameter estimates
            if save_results == 1

                if model_type == 1
                    save(['SF_compounds_', num2str(mod_name), '_sim_data.mat'], 'sim_data')% save

                elseif model_type == 2

                    save(['SF_features_', num2str(mod_name), '_sim_data.mat'], 'sim_data')% save

                elseif model_type == 5

                    save(['SF_features_transfer_choice_mixing_', num2str(mod_name), '_sim_data.mat'], 'sim_data')% save

                end
            end

            cd ..
          
        end
    end
end
