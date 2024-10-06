function fit_models_main(subj, iter, model_type, multi_lr, multi_tau, prior_sigma, MAP_estimation, oracle, run_cluster)

%% Prepare data sets and fit different Successor Feature Models (using Maximum Likelihood estimation, fmincon)
% Models: state-state/compound-compound, feature-feature, or hybrid model to data
% Models feature representation (average expected, discounted future
% feature/state occupancy) and finds optimal parameters (learning rate and
% softmax stochasticity parameter, updating parameter of weighted average between representations and 
% initial bias towards one or the other representation at first trial)

% clear all; close all; clc
disp(['Fitting model: ' num2str(model_type)])
disp(['Subject: ' num2str(subj)])

format compact

script_dir = pwd;

if run_cluster == 1
    source_path = '/data/holly-host/lluettgau/review_code_data';
else
    source_path = '/Users/lenna/Desktop/Work/PostDoc/Projects/Markov/task/review_code_data';
end

%% settings
rng((subj * iter)) %set rng seed

fitted_initial_values = 1; %use MLEs from SF feature model as starting points for SF feature transfer model
fixed_parameters = 0;

if MAP_estimation == 1
    opt_name = 'MAP';
else
    opt_name = 'MLE';
end

if model_type == 5
    if fitted_initial_values == 1
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

        all_params = load([source_path, filesep, 'fitted_parameters', filesep, 'fitted_parameters_', opt_name, '_SF_features_', num2str(mod_name), '_rand_inival.mat']); %specify folder for fitted params here
        ini_all = all_params.comp_params.params;
        num_choices_tmp = all_params.comp_params.nsamples;
        subIDs_all = all_params.comp_params.subIDs;
        ini_all = [subIDs_all ini_all];
    end
end

%temporal discounting parameter
gamma = .9;

%% specify data and paths 
res_dir = [source_path, filesep, 'fitted_parameters', filesep, 'cluster']; %specify folder for fitted params here
data_dir = [source_path, filesep, 'data']; %specify data folder here
addpath(genpath([source_path, filesep, 'functions']))

cd(data_dir)

task_data = readtable('all_data.csv');
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

%% loop to prepare data and fit subjects
cd(script_dir)

clear cp 

feature_names = {'A', 'B', 'C', 'D', ...
                 '1', '2', '3', '4', '5', '6'}';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%% Preprocessing/preparation of data
preproc_outputs = preprocess_data(task_data, task_begin_idx, task_end_idx, subj, feature_names, model_type);
condition = str2double(preproc_outputs.experimental_data.Spreadsheet{1}(21));
subID = preproc_outputs.experimental_data.ParticipantPrivateID(1);

%prepare inputs to model
preproc_outputs.condition = condition;
preproc_outputs.fixed_params = fixed_parameters;
preproc_outputs.gamma = gamma;
preproc_outputs.simulation = 0;
preproc_outputs.multi_lr = multi_lr;
preproc_outputs.multi_tau = multi_tau;

preproc_outputs.oracle = oracle;  


%% Model fitting
%% SF compound model
if model_type == 1
    
    %function optimization
    x = zeros(1,multi_lr+multi_tau);
    lb = zeros(1,numel(x));
    if multi_lr == 1
        if multi_tau == 1
            ub = [1, 100];
            %for MAP, prior simulation parameters
            mu = [0.1; .3];
        else
            ub = [1, 100, 100];
            %for MAP, prior simulation parameters
            mu = [0.1; .3; .6];
        end
    elseif multi_lr == 2
        if multi_tau == 1
            ub = [1, 1, 100];
            %for MAP, prior simulation parameters
            mu = [0.01; 0.1; .3];
        else
            ub = [1, 1, 100, 100];
            %for MAP, prior simulation parameters
            mu = [0.01; 0.1; .3; .6];
        end

    end


    %set up MAP estimation
    nui = eye(numel(x)) * prior_sigma;
    np = length(mu);
    log_prior = @(x) -1/2 * (x - mu)' * nui * (x - mu) - np/2 * log(2*pi) - 1/2 * log(1/det(nui));

    %define functions to be optimized
    fun_MAP = @(x) SF_compound_model_all(preproc_outputs, x) - log_prior(x'); %add prior penalty for MAP
    fun_MLE = @(x) SF_compound_model_all(preproc_outputs, x); %MLE function

    options = optimoptions('fmincon','Display','off');
    
    %repeat iter times for different initial values
    ini = lb + (ub - lb) .* rand(size(lb));

    if MAP_estimation == 1
        [minimum,energy_map,EXITFLAG,~,~,~,hessian] = fmincon(fun_MAP, ini,[],[],[],[],lb,ub,[],options);
    else
        [minimum,energy,EXITFLAG,~,~,~,hessian] = fmincon(fun_MLE, ini,[],[],[],[],lb,ub,[],options);
    end

    %use optimal MAP parameters to obtain LLE
    %get number of choices
    [energy, model_outputs] = SF_compound_model_all(preproc_outputs, minimum);
    num_all_choices = model_outputs.num_all_choices;

 %% SF feature model
elseif model_type == 2

    %prepare inputs to model
    preproc_outputs.transfer_model = 0;

    %function optimization
    x = zeros(1,multi_lr+multi_tau);
    lb = zeros(1,numel(x));

    if multi_lr == 1
        if multi_tau == 1
            ub = [1, 100];
            %for MAP, prior simulation parameters
            mu = [0.1; .3];
        else
            ub = [1, 100, 100];
            %for MAP, prior simulation parameters
            mu = [0.1; .3; .6];
        end
    elseif multi_lr == 2
        if multi_tau == 1
            ub = [1, 1, 100];
            %for MAP, prior simulation parameters
            mu = [0.01; 0.1; .3];
        else
            ub = [1, 1, 100, 100];
            %for MAP, prior simulation parameters
            mu = [0.01; 0.1; .3; .6];
        end

    end

    %set up MAP estimation
    nui = eye(numel(x)) * prior_sigma;
    np = length(mu);
    log_prior = @(x) -1/2 * (x - mu)' * nui * (x - mu) - np/2 * log(2*pi) - 1/2 * log(1/det(nui));

    %define functions to be optimized
    fun_MAP = @(x) SF_feature_model_all(preproc_outputs, x) - log_prior(x'); %add prior penalty for MAP
    fun_MLE = @(x) SF_feature_model_all(preproc_outputs, x); %MLE function

    options = optimoptions('fmincon','Display','off');
    
    %repeat iter for different initial values
    ini = lb + (ub - lb) .* rand(size(lb));

    if MAP_estimation == 1

        [minimum,energy_map,EXITFLAG,~,~,~,hessian] = fmincon(fun_MAP, ini,[],[],[],[],lb,ub,[],options);

    else
        [minimum,energy,EXITFLAG,~,~,~,hessian] = fmincon(fun_MLE, ini,[],[],[],[],lb,ub,[],options);
    end

    %use optimal MAP parameters to obtain LLE
    %get number of choices
    [energy, model_outputs] = SF_feature_model_all(preproc_outputs, minimum);
    num_all_choices = model_outputs.num_all_choices;

    %% SF feature model with transfer learning, separate learning rate and mixing parameter
elseif model_type == 5

    %prepare inputs to model
    preproc_outputs.transfer_model = 1;

    %function optimization
    x = zeros(1,multi_lr + multi_tau + 1);
    lb = zeros(1,numel(x));

    if multi_lr == 1
        if multi_tau == 1
            ub = [1, 100, 1];
            %for MAP, prior simulation parameters
            mu = [0.1; .3; .25];
        else
            ub = [1, 100, 100, 1];
            %for MAP, prior simulation parameters
            mu = [0.1; .3; .6; .25];
        end
    elseif multi_lr == 2
        if multi_tau == 1
            ub = [1, 1, 100, 1];
            %for MAP, prior simulation parameters
            mu = [0.01; 0.1; .3; .25];
        else
            ub = [1, 1, 100, 100, 1];
            %for MAP, prior simulation parameters
            mu = [0.01; 0.1; .3; .6; .25];
        end

    end

    %set up MAP estimation
    nui = eye(numel(x)) * prior_sigma;
    np = length(mu);
    log_prior = @(x) -1/2 * (x - mu)' * nui * (x - mu) - np/2 * log(2*pi) - 1/2 * log(1/det(nui));

    %define functions to be optimized
    fun_MAP = @(x) SF_feature_model_all(preproc_outputs, x) - log_prior(x'); %add prior penalty for MAP
    fun_MLE = @(x) SF_feature_model_all(preproc_outputs, x); %MLE function

    options = optimoptions('fmincon','Display','off');
        
    ini_idx = find(ini_all(:,1) == subID);

    %repeat iter times for different initial values
    if fitted_initial_values == 0
        ini = lb + (ub - lb) .* rand(size(lb));
    else
        if multi_lr == 1
            if multi_tau == 1
                ini_w = lb(3) + (ub(3) - lb(3)) .* rand(size(lb(3)));
            else
                ini_w = lb(4) + (ub(4) - lb(4)) .* rand(size(lb(4)));
            end
        else
            if multi_tau == 1
                ini_w = lb(4) + (ub(4) - lb(4)) .* rand(size(lb(4)));
            else
                ini_w = lb(5) + (ub(5) - lb(5)) .* rand(size(lb(5)));
            end
        end
        
        ini = [ini_all(ini_idx,2:end) ini_w];

    end

    if MAP_estimation == 1

        [minimum,energy_map,EXITFLAG,~,~,~,hessian] = fmincon(fun_MAP, ini,[],[],[],[],lb,ub,[],options);

        %use optimal MAP parameters to obtain LLE
        [energy, model_outputs] = SF_feature_model_all(preproc_outputs, minimum);
    else
        [minimum,energy,EXITFLAG,~,~,~,hessian] = fmincon(fun_MLE, ini,[],[],[],[],lb,ub,[],options);

    end
    %get number of choices
    num_all_choices = num_choices_tmp(ini_idx,1);

end %end model_type statement

% Save results for this subject and iteration
results = [ini, subID, subj, iter, minimum, energy, exp(-energy/num_all_choices), num_all_choices];

if oracle == 0
    %save parameter estimates
    save([res_dir, filesep, 'fitted_params_' opt_name '_model_' num2str(model_type) ...
        '_alpha_' num2str(multi_lr) '_tau_' num2str(multi_tau) '_subject_' num2str(subj) '_iteration_' num2str(iter) '.mat'],'results')
else
    %save parameter estimates
    save([res_dir, filesep, 'fitted_params_' opt_name '_oracle_model_' num2str(model_type) ...
        '_alpha_' num2str(multi_lr) '_tau_' num2str(multi_tau) '_subject_' num2str(subj) '_iteration_' num2str(iter) '.mat'],'results')
end

end