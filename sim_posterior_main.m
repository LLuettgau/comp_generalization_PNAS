%% Posterior simulations: Successor Compound/Feature Models - 
% use fitted parameter values to generate choices

clear all; close all; clc

format compact

addpath(genpath(['.', filesep, 'functions']))
addpath(genpath(['.', filesep, 'posterior_simulations']))
addpath(genpath('.'))

%% settings
% model_type = 1; %Compounds (1), 
% % Features (2), 
% % Features (transfer using permutations, choice mixing param, 5)
% multi_lr = 1; %1 or 2 learning rates (prior/transfer separately)

save_data = 1;
save_plots = 0;
fixed_parameters = 0;
oracle = 0;
MAP_estimation = 0;

if MAP_estimation == 1
    if oracle == 1
        opt_name = 'MAP_oracle';
    else
        opt_name = 'MAP';
    end
else
    if oracle == 1
        opt_name = 'MLE_oracle';
    else
        opt_name = 'MLE';
    end
end

niter = 1000; %how many simulations, 1000 in main analysis

%start parallel processing
delete(gcp('nocreate'))

num_workers = 8;
pool = parpool(num_workers);

%% specify data, sample and paths
res_dir = ['.', filesep, 'posterior_simulations']; %specify folder for fitted params here
data_dir = ['.', filesep, 'data']; %specify data folder here
plotpath = ['.', filesep, 'posterior_simulations', filesep, 'plots']; %specify folder for plots here
modeldir = ['.', filesep, 'fitted_parameters/'];

%get task data
task_data = load_data(data_dir);

%remove unnecessary variables
columns_to_remove = {'background', 'Pause_duration', 'ITI_duration', 'nextRow', ...
    'event', 'randomise_blocks','YCoordinate','XCoordinate', ...
    'Dishonest', 'Incorrect','Attempt','ZoneType', ...
    'ZoneName', 'ScreenNumber','Attempt','SpreadsheetRow', ...
    'PercentageScoreTransfer', 'IncorrectScoreTransfer','CorrectScoreTransfer','PercentageScore', ...
    'IncorrectScore', 'CorrectScore','checkpoint_k8iu','checkpoint_v7br', ...
    'checkpoint_mgdx', 'branch_rctx','randomiser_xx15','TaskVersion', ...
    'TaskName', 'Checkpoint','ParticipantViewportSize','ParticipantMonitorSize', ...
    'ParticipantBrowser', 'ParticipantOS','ParticipantDevice','ParticipantDeviceType', ...
    'ParticipantCompletionCode', 'ParticipantStatus','ParticipantStartingGroup','ParticipantPublicID', ...
    'ScheduleID', 'RepeatKey','TreeNodeKey','ExperimentVersion', ...
    'ExperimentID', 'LocalDate','LocalTimezone','LocalTimestamp', ...
    'UTCDate', 'UTCTimestamp'};

task_data = removevars(task_data, columns_to_remove);

task_begin_idx = find(strcmp(task_data.display, 'instructionsExperiment'));
task_end_idx = find(strcmp(task_data.display, 'End'));

cd ..

feature_names = {'A', 'B', 'C', 'D', ...
    '1', '2', '3', '4', '5', '6'}';


%get subIDs
subIDs_models_tmp = load([modeldir, 'subIDs_cell.mat']);
subIDs_models = subIDs_models_tmp.subject_ID_cell{1,1}(:);
subIDs_models = sort(subIDs_models);

%specify gamma
gamma = .9;

for model_type = 5%[1 2 5]
    for multi_lr = 1:2
        for multi_tau = 2%1:2

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
            
            model_name = ['fitted_parameters_', opt_name, '_SF_compounds_', num2str(mod_name), '_rand_inival.mat'];

            transfer_model = 0;

        elseif model_type == 2
            
            model_name = ['fitted_parameters_', opt_name, '_SF_features_', num2str(mod_name), '_rand_inival.mat'];

            transfer_model = 0;

        elseif model_type == 5
 
            model_name = ['fitted_parameters_', opt_name, '_SF_features_transfer_', num2str(mod_name), '_1choiceMix_fitted_inival.mat'];
            
            transfer_model = 1;
        end

        ind_data_points = 0; %plot individual data points?


        %% get fitted parameters
        params_tmp = importdata([modeldir,model_name]);
        fitted_params = params_tmp.params;
        fitted_params = [subIDs_models fitted_params];


        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Model simulation

        tic
        for subj = 1:length(task_begin_idx)

            %% Preprocessing/preparation of data
            preproc_outputs = preprocess_data(task_data, task_begin_idx, task_end_idx, subj, feature_names, model_type);

            %read experimental data to get schedule info
            experimental_data = preproc_outputs.experimental_data;
            schedule_name = experimental_data.Spreadsheet(1);

            subject_ID(subj,1) = experimental_data.ParticipantPrivateID(1);
            condition(subj,1) = str2double(schedule_name{1}(21));

            param_idx = find(fitted_params(:,1) == subject_ID(subj,1));

            %get fitted params
            params = fitted_params(param_idx,2:size(fitted_params,2));

            %prepare inputs to model
            preproc_outputs.fixed_params = fixed_parameters;
            preproc_outputs.gamma = gamma;
            preproc_outputs.condition = condition(subj);
            preproc_outputs.transfer_model = transfer_model;
            preproc_outputs.simulation = 1;
            preproc_outputs.multi_lr = multi_lr;
            preproc_outputs.multi_tau = multi_tau;
            preproc_outputs.oracle = oracle;


            %% simualtions using fit params
            parfor iter = 1:niter

                % preproc_outputs = preprocess_data(task_data, task_begin_idx, task_end_idx, subj, feature_names, model_type);

                local_sim_results = nan(1, 10);

                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %SF compound model
                if model_type == 1

                    %run simulation
                    [~, model_outputs] = SF_compound_model_all(preproc_outputs, params);

                    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %SF feature model
                elseif model_type == 2

                    %run simulation
                    [~, model_outputs] = SF_feature_model_all(preproc_outputs, params);

                    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %SF feature model (transfer)
                elseif model_type == 5

                    %run simulation
                    [~, model_outputs] = SF_feature_model_all(preproc_outputs, params);

                end


                local_sim_results(1) = subject_ID(subj, 1);
                local_sim_results(2) = condition(subj, 1);
                %exp probes, first factor, prior learning
                local_sim_results(3) = nanmean(model_outputs.prior_sim_res_corr_exp(preproc_outputs.exp_trial_type==1));
                %exp probes, second factor, prior learning
                local_sim_results(4) = nanmean(model_outputs.prior_sim_res_corr_exp(preproc_outputs.exp_trial_type==2));
                %inf probes, first factor, prior learning
                local_sim_results(5) = nanmean(model_outputs.prior_sim_res_corr_inf(preproc_outputs.inf_trial_type==1));
                %inf probes, second factor, prior learning
                local_sim_results(6) = nanmean(model_outputs.prior_sim_res_corr_inf(preproc_outputs.inf_trial_type==2));

                %exp probes, first factor, transfer learning
                local_sim_results(7) = nanmean(model_outputs.trans_sim_res_corr_exp(preproc_outputs.trans_exp_trial_type==1));
                %exp probes, second factor, transfer learning
                local_sim_results(8) = nanmean(model_outputs.trans_sim_res_corr_exp(preproc_outputs.trans_exp_trial_type==2));
                %inf probes, first factor, transfer learning
                local_sim_results(9) = nanmean(model_outputs.trans_sim_res_corr_inf(preproc_outputs.trans_inf_trial_type==1));
                %inf probes, second factor, transfer learning
                local_sim_results(10) = nanmean(model_outputs.trans_sim_res_corr_inf(preproc_outputs.trans_inf_trial_type==2));

                %store results
                intermediate_results{iter} = local_sim_results;

            end %end of iteration loop

            % Store subject_results in the all_subject_results cell array
            all_subject_results{subj} = intermediate_results;

        end %end subject loop


        % Combine the intermediate results into sim_results after the parfor loop
        for subj = 1:length(task_begin_idx)
            for iter = 1:niter
                sim_results(subj, :, iter) = all_subject_results{subj}{iter};
            end
        end

        toc

        if model_type == 5 && multi_lr == 2
            delete(pool)
        end


        %% Plots

        %averages in different trial types
        res_to_plot = [(nanmean(sim_results(:,3,:),3) + nanmean(sim_results(:,4,:),3))/2 ...
            (nanmean(sim_results(:,5,:),3) + nanmean(sim_results(:,6,:),3))/2 ...
            (nanmean(sim_results(:,7,:),3) + nanmean(sim_results(:,8,:),3))/2 ...
            (nanmean(sim_results(:,9,:),3) + nanmean(sim_results(:,10,:),3))/2];

        %experience
        exp_res_to_plot = [(nanmean(sim_results(:,3,:),3) + nanmean(sim_results(:,4,:),3))/2 ...
            nanmean(sim_results(:,3,:),3) nanmean(sim_results(:,4,:),3) ...
            (nanmean(sim_results(:,7,:),3) + nanmean(sim_results(:,8,:),3))/2 ...
            nanmean(sim_results(:,7,:),3) nanmean(sim_results(:,8,:),3)];


        %inference
        inf_to_plot = [(nanmean(sim_results(:,5,:),3) + nanmean(sim_results(:,6,:),3))/2 ...
            nanmean(sim_results(:,5,:),3) nanmean(sim_results(:,6,:),3) ...
            (nanmean(sim_results(:,9,:),3) + nanmean(sim_results(:,10,:),3))/2 ...
            nanmean(sim_results(:,9,:),3) nanmean(sim_results(:,10,:),3)];



        diffs(:,1) = nanmean(nanmean(sim_results(condition==1,7,:),3)) - nanmean(nanmean(sim_results(condition==2,7,:),3));
        diffs(:,2) = nanmean(nanmean(sim_results(condition==1,8,:),3)) - nanmean(nanmean(sim_results(condition==2,8,:),3));
        diffs(:,3) = nanmean(nanmean(sim_results(condition==1,9,:),3)) - nanmean(nanmean(sim_results(condition==2,9,:),3));
        diffs(:,4) = nanmean(nanmean(sim_results(condition==1,10,:),3)) - nanmean(nanmean(sim_results(condition==2,10,:),3));

        for k = 1:4
            diffs(2,k) = nanstd(nanmean(sim_results(condition==1,6+k,:),1)) / sqrt(length(nanmean(sim_results(condition==1,6+k,:),1))) + ...
                nanstd(nanmean(sim_results(condition==2,6+k,:),1)) / sqrt(length(nanmean(sim_results(condition==2,6+k,:),1)));
        end

        %% Plotting
        ylim_lower = 0.4;
        y_ticks = ylim_lower:.2:1;
        results(:,8) = condition;
        results(:,1) = sim_results(:,1,1);

        plotting_simulations

        if save_plots == 1
            screen_size = get(0, 'ScreenSize');
            origSize = get(h, 'Position'); % grab original on screen size
            set(h, 'Position', [0 0 screen_size(3) screen_size(4)/2 ] ); %set to scren size
            set(h,'PaperPositionMode','auto') %set paper pos for printing
            % saveas(h, 'fig_trajectories.png') % save figure
            ax = h;
            exportgraphics(ax,[plotpath, filesep, 'fig_pcorrect_all_' model_name(19:end-4) '.png'],'Resolution',300)
            set(h,'Position', origSize) %set back to original dimensions
        end


        %% condition differences of interest

        positions = [1 2 3 4 5];
        pos = [positions-.3; ...
            positions+.3];
        size_vec = [1 2 3 4 5];


        h = figure('visible','on');
        b = bar(1:4,diffs(1,:),'FaceColor','flat', 'BarWidth', 0.6, 'linewidth', 2.5, 'FaceAlpha',.5); hold all;
        b.CData(1,:) = cols.green;
        b.CData(2,:) = cols.k;
        b.CData(3,:) = cols.green;
        b.CData(4,:) = cols.k;

        errorbar(positions(1),nanmean(diffs(1,1)), diffs(2,1), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
        errorbar(positions(2),nanmean(diffs(1,2)), diffs(2,2), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
        errorbar(positions(3),nanmean(diffs(1,3)), diffs(2,3), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
        errorbar(positions(4),nanmean(diffs(1,4)), diffs(2,4), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);

        set(findobj(gca,'type','line'),'linew',2)

        %ylim([0 1]);
        ylim([-0.1 .1]);
        %xlim([0 5]);
        xlim([0.5 4.5]);
        box off
        set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
        ybounds = ylim;
        yline(0,'--','linewidth',2.5)
        %     set(gca,'YTick',[0 .25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
        set(gca,'YTick',-0.1:.05:.1, 'FontSize',26,'FontName', 'Arial');
        set(gca,'TickDir','out')
        set(gca,'xtick',1:4)

        set(gcf,'color','w');
        set(gca,'ycolor',cols.k)
        set(gca,'xcolor',cols.k)

        % prepend a color for each tick label
        LabelsY = get(gca,'YTickLabel');
        ticklabels_ynew = cell(size(LabelsY));
        for i = 1:length(LabelsY)
            ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
        end
        % set the tick labels
        set(gca, 'YTickLabel', ticklabels_ynew,'FontSize',20,'FontName', 'Arial');
        ylabel('\Delta Probability Correct','FontSize',20,'FontName', 'Arial')

        set(gca,'XTickLabel', Labels_trans, 'FontSize',26,'FontName', 'Arial');
        set(gca,'XTickLabelRotation', 45);
        ticklabels_new = cell(size(Labels_trans));
        for i = 1:length(Labels_trans)
            ticklabels_new{i} = ['\color{black} ' Labels_trans{i}];
        end
        % set the tick labels
        set(gca, 'XTickLabel', ticklabels_new,'FontSize',26,'FontName', 'Arial');


        %% heatmap of performance, per condition
        h = figure

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
            exportgraphics(ax,['heatmap_posterior_simulation_model_' num2str(model_type) '.png'],'Resolution',600)
            set(h,'Position', origSize) %set back to original dimensions
        end

        cd ..

        %% save data for plotting
        sim_data.results = results;
        sim_data.res_to_plot = res_to_plot;
        sim_data.exp_res_to_plot = exp_res_to_plot;
        sim_data.inf_to_plot = inf_to_plot;
        sim_data.diffs = diffs;

        if save_data == 1
            if model_type == 1

                save(['SF_compounds_', opt_name, '_', num2str(mod_name), '_ppc_sim_data.mat'], 'sim_data')% save

            elseif model_type == 2

                save(['SF_features_', opt_name, '_', num2str(mod_name), '_ppc_sim_data.mat'], 'sim_data')% save

            elseif model_type == 5
                save(['SF_features_transfer_', opt_name, '_', num2str(mod_name), '_ppc_sim_data.mat'], 'sim_data')% save

            end
        end

        cd ..
        end
    end
end
