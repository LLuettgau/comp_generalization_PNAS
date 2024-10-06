restoredefaultpath
addpath(genpath(['.', filesep, 'posterior_simulations']))
addpath(genpath(['.', filesep, 'functions']))
addpath(genpath('.'))

res_dir = ['.', filesep, 'posterior_simulations']; %specify folder for fitted params here
plotpath = ['.', filesep, 'posterior_simulations/plots/']; 
%specify colors for plots
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
cols.purple = [151 99 179]/255;
cols.dpurple = [0.400 0.100 0.500];
cols.dyellow = [0.9290 0.6940 0.1250];
cols.dred = [0.6350 0.0780 0.1840];
cols.lgreen = [0.4660 0.6740 0.1880];

cols.dblue = [60 54 152]/255;
cols.magenta = [178 45 110]/255;
cols.dorange = [240 86 32]/255;
cols.gold = [252, 183, 30]/255;

cols.cream = [220 190 110]/255;
cols.lyellow = [254 153 41]/255;
cols.medorange = [217 95 14]/255;
cols.dred = [153 52 4]/255;


ylim_lower = 0.4;
y_ticks = ylim_lower:.2:1;

Labels_trans = {'4-cycle prior - 4-cycle', '6-cycle prior - 4-cycle', '4-cycle prior - 6-cycle', '6-cycle prior - 6-cycle'};
Labels_diff_trans = {'Experience: 4-cycle', ...
                     'Experience: 6-cycle', ...
                     'Inference: 4-cycle', ...
                     'Inference: 6-cycle'};
Labels_within_diff_trans = {'Experience: 4-cycle prior', ...
                            'Experience: 6-cycle prior', ...
                            'Inference: 4-cycle prior', ...
                            'Inference: 6-cycle prior'};



sim_file_names = {'SF_compounds_MLE_1Alpha_1Tau_ppc_sim_data.mat', ...
                  'SF_compounds_MLE_1Alpha_2Tau_ppc_sim_data.mat', ...
                  'SF_compounds_MLE_2Alpha_1Tau_ppc_sim_data.mat', ...
                  'SF_compounds_MLE_2Alpha_2Tau_ppc_sim_data.mat', ...
                  ...
                  'SF_features_MLE_1Alpha_1Tau_ppc_sim_data.mat', ...
                  'SF_features_MLE_1Alpha_2Tau_ppc_sim_data.mat', ...
                  'SF_features_MLE_2Alpha_1Tau_ppc_sim_data.mat', ...
                  'SF_features_MLE_2Alpha_2Tau_ppc_sim_data.mat', ...
                  ...
                  'SF_features_transfer_MLE_1Alpha_1Tau_ppc_sim_data.mat', ...
                  'SF_features_transfer_MLE_1Alpha_2Tau_ppc_sim_data.mat', ...
                  'SF_features_transfer_MLE_2Alpha_1Tau_ppc_sim_data.mat', ...
                  'SF_features_transfer_MLE_2Alpha_2Tau_ppc_sim_data.mat'};

model_names = {'1\alpha, 1\tau', ...
               '1\alpha, 2\tau', ...
               '2\alpha, 1\tau', ...
               '2\alpha, 2\tau'};

%% Observed data
mod_num = 1;

load([res_dir, filesep, 'observed_data.mat'])
results = obs_data.results;
res_to_plot = obs_data.res_to_plot;
exp_res_to_plot = obs_data.sl_res_to_plot;
inf_to_plot = obs_data.inf_to_plot;

results_sorted = sortrows(results,1);
subIDs = results_sorted(:,8);


h = figure('visible','on');
t = tiledlayout(4,4,'TileSpacing','loose', 'Padding', 'tight');
counter = 1;

for i = 1:4
    all_sim_plots
end


%% plot model predictions
mod_num = 3;

for k = 1:size(sim_file_names,2)
    %load data
    clear all_sim_data
    load([res_dir, filesep, sim_file_names{k}]);

    results = sim_data.results;
    res_to_plot = sim_data.res_to_plot;
    exp_res_to_plot = sim_data.exp_res_to_plot;
    inf_to_plot = sim_data.inf_to_plot;

    %make plots
    all_sim_plots
    if sum(k == 1:4) > 0 
        title(model_names{k},'FontSize',24,'FontName', 'Arial');
    end

    if sum(k == 9:12) > 0
        
        set(gca,'XTickLabel', Labels_diff_trans, 'FontSize',20,'FontName', 'Arial');
        set(gca,'XTickLabelRotation', 45);
        ticklabels_new = cell(size(Labels_diff_trans));
        for i = 1:length(Labels_diff_trans)
            ticklabels_new{i} = ['\color{black} ' Labels_diff_trans{i}];
        end
        % set the tick labels
        set(gca, 'XTickLabel', ticklabels_new,'FontSize',20,'FontName', 'Arial');
    
    end

    plotting_simulations
end
    
% Create shared ylabel
ylabel(t,'\Delta Probability Correct', 'FontSize',24,'FontName', 'Arial')





% h = figure('visible','on');
% t = tiledlayout(1,4,'TileSpacing','loose', 'Padding', 'tight');
% counter = 1;
% 
% 
% 
% Labels_prior = {'4-cycle prior - 4-cycle',  '6-cycle prior - 4-line', '4-cycle prior - 6-line', '6-cycle prior - 6-cycle'};
% Labels_trans = {'4-cycle prior - 4-cycle', '6-cycle prior - 4-cycle', '4-cycle prior - 6-cycle', '6-cycle prior - 6-cycle'};
% 
% 
% for k = 10%1:size(sim_file_names,2)
%     %load data
%     clear all_sim_data
%     load([res_dir, filesep, sim_file_names{k}]);
% 
%     results = sim_data.results;
%     res_to_plot = sim_data.res_to_plot;
%     exp_res_to_plot = sim_data.exp_res_to_plot;
%     inf_to_plot = sim_data.inf_to_plot;
% 
%     nexttile
%     %prior learning - experience probes
%     avg_bar_scatter_simulations(results, exp_res_to_plot, Labels_prior, cols, [2 3], 1, ylim_lower, 4)
% 
%     %prior learning - inference probes
%     nexttile
%     avg_bar_scatter_simulations(results, inf_to_plot, Labels_prior, cols, [2 3], 1, ylim_lower, 4)
% 
%     %transfer learning - experience probes
%     nexttile
%     avg_bar_scatter_simulations(results, exp_res_to_plot, Labels_trans, cols, [5 6], 2, ylim_lower, 4)
% 
%     %transfer learning - inference probes
%     nexttile
%     avg_bar_scatter_simulations(results, inf_to_plot, Labels_trans, cols, [5 6], 2, ylim_lower, 4)
% 
% end
% %% Structural inference model
% mod_num = 2;
% 
% all_sim_plots


% screen_size = get(0, 'ScreenSize');
% origSize = get(h, 'Position'); % grab original on screen size
% set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size
% set(h,'PaperPositionMode','auto') %set paper pos for printing
% % saveas(h, 'fig_trajectories.png') % save figure
exportgraphics(t,[plotpath 'posterior_simulations_all.png'],'Resolution',1200)
set(h,'Position', origSize) %set back to original dimensions

