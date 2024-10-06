%% Prepare data sets and fit different Successor Feature Models (using Maximum Likelihood estimation, fmincon)
% Models: state-state/compound-compound, feature-feature, or hybrid model to data
% Models feature representation (average expected, discounted future
% feature/state occupancy) and finds optimal parameters (learning rate and
% softmax stochasticity parameter, updating parameter of weighted average between representations and 
% initial bias towards one or the other representation at first trial)

clear all; close all; clc

format compact

script_dir = pwd;


%% settings
% model_type = 1; %(1) SF compound model (compound-compund)
%                 %(2) SF feature model (feature-feature)
%                 %(3) SF hybrid (weighted average between SF compound and SF feature model)
%                 %(5) SF feature model (feature-feature) with transfer learning , separate learning rate and mixing parameter

MAP_estimation = 0;

for oracle = 1

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

    for model_type = 5%[1 2 5]
        for multi_lr = 2
            for multi_tau = 2

                clear cp_tmp
                clear cp

                %% specify data and paths
                res_dir = ['.', filesep, 'fitted_parameters']; %specify folder for fitted params here

                if multi_tau == 1 && model_type < 5
                    all_params_dir = dir(['.', filesep, 'fitted_parameters', filesep, 'cluster', filesep, ...
                        'fitted_params_', opt_name, '_model_', num2str(model_type), ...
                        '_alpha_', num2str(multi_lr), '_subject*_iteration*.mat']); %specify folder for fitted params here
                elseif multi_tau == 2 || model_type == 5
                    all_params_dir = dir(['.', filesep, 'fitted_parameters', filesep, 'cluster', filesep, ...
                        'fitted_params_', opt_name, '_model_', num2str(model_type), ...
                        '_alpha_', num2str(multi_lr), '_tau_', num2str(multi_tau), '*_iteration*.mat']); %specify folder for fitted params here
                end

                %% restructure data

                % Loop through the files and extract subject numbers
                for i = 1:length(all_params_dir)
                    % Get the filename
                    filename = all_params_dir(i).name;

                    % Extract the subject number from the filename
                    sub_iter = regexp(filename, '_subject_(\d+)_iteration_(\d+).mat', 'tokens');

                    loadedData = load([all_params_dir(i).folder, filesep, filename]);
                    loadedData = loadedData.results;

                    cp_tmp.params(i,1:2) = str2double(sub_iter{1});
                    cp_tmp.params(i,3:size(loadedData,2)+2) = loadedData;
                end

                %remove initial values
                if model_type < 5
                    cp_tmp.params(:,3:3 + multi_lr + multi_tau - 1) = [];
                else
                    cp_tmp.params(:,3:4 + multi_lr + multi_tau - 1) = [];
                end

                subIDs = cp_tmp.params(:,3);
                cp_tmp.params(:,3) = [];
                cp_tmp.params = [subIDs cp_tmp.params];

                %ascertain match between name and save sub/iter
                if sum(cp_tmp.params(:,2) - cp_tmp.params(:,4) == 0) == length(cp_tmp.params)
                    % Extract unique subject IDs from the first column
                    sub_nums = unique(cp_tmp.params(:, 1));

                    % Initialize a matrix to store results
                    cp.params = zeros(length(sub_nums), size(cp_tmp.params, 2));

                    % Iterate over unique subject IDs
                    for i = 1:length(sub_nums)
                        % Extract data for the current subject
                        subjectData = cp_tmp.params(cp_tmp.params(:, 1) == sub_nums(i), :);

                        % Find the row index with the minimum value in LLE
                        % column
                        [~, minRowIndex] = min(subjectData(:, end-2));

                        % Store the row with the minimum value in resultMatrix
                        cp.params(i, :) = subjectData(minRowIndex, :);
                    end

                else
                    disp('no match between unique subject IDs and file names!')
                end

                cp.params(:,2:5) = [];
                cp.params = sortrows(cp.params,1);
                subIDs = cp.params(:,1);

                cp.params(:,1) = [];

                nanmedian(cp.params)
                nanmean(cp.params)
                min(cp.params)
                max(cp.params)
                % 
                % figure; hist(cp.params(:,1))
                % figure; hist(cp.params(:,2))
                % figure; hist(cp.params(:,3))
                % figure; hist(cp.params(:,4))
                % figure; hist(cp.params(:,5))

                %write into struct
                comp_params.LLE = cp.params(:,end-2);
                comp_params.ntrials_explained = cp.params(:,end-1);
                comp_params.nsamples = cp.params(:,end);
                comp_params.params = cp.params(:,1:end-3);
                comp_params.subIDs = subIDs;

                %save parameter estimates
                cd(res_dir)
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

                    save([pwd,filesep,'fitted_parameters_', opt_name, '_SF_compounds_', num2str(mod_name), '_rand_inival.mat'],'comp_params');

                elseif model_type == 2

                    save([pwd,filesep,'fitted_parameters_', opt_name, '_SF_features_', num2str(mod_name), '_rand_inival.mat'],'comp_params');

                elseif model_type == 5

                    save([pwd,filesep,'fitted_parameters_', opt_name, '_SF_features_transfer_', num2str(mod_name), '_1choiceMix_fitted_inival.mat'],'comp_params');
                end

                %%% move used parameters to processed params folder
                % cd cluster/
                % sourceDirectory = pwd;  % specify the source directory
                % destinationDirectory = [pwd,filesep,'processed_params/'];  % specify the destination directory
                % 
                % for i = 1:length(all_params_dir)
                %     filename = all_params_dir(i).name;
                %     bashCommand = ['mv ' sourceDirectory '/' filename ' ' destinationDirectory];
                %     system(bashCommand);
                % end

                cd ..
            end
        end
    end
end