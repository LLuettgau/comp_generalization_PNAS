restoredefaultpath
addpath(genpath(['.', filesep, 'prior_simulations']))
addpath(genpath(['.', filesep, 'functions']))
addpath(genpath('.'))

res_dir = ['.', filesep, 'prior_simulations']; %specify folder for fitted params here
plotpath = ['.', filesep, 'prior_simulations/plots/']; %specify folder for fitted params here

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

%% combine plots across simulations
% sim_file_names = {'SF_compounds_1Alpha_1Tau_sim_data.mat', ...
%                   'SF_compounds_1Alpha_2Tau_sim_data.mat', ...
%                   'SF_compounds_2Alpha_1Tau_sim_data.mat', ...
%                   'SF_compounds_2Alpha_2Tau_sim_data.mat', ...
%                   ...
%                   'SF_features_1Alpha_1Tau_sim_data.mat', ...
%                   'SF_features_1Alpha_2Tau_sim_data.mat', ...
%                   'SF_features_2Alpha_1Tau_sim_data.mat', ...
%                   'SF_features_2Alpha_2Tau_sim_data.mat', ...
%                    ...
%                   'SF_features_transfer_choice_mixing_1Alpha_1Tau_sim_data.mat', ...
%                   'SF_features_transfer_choice_mixing_1Alpha_2Tau_sim_data.mat', ...
%                   'SF_features_transfer_choice_mixing_2Alpha_1Tau_sim_data.mat', ...
%                   'SF_features_transfer_choice_mixing_2Alpha_2Tau_sim_data.mat'};

% model_names = {'SF_C 1\alpha, 1\tau', ...
%                  'SF_C 1\alpha, 2\tau', ...
%                  'SF_C 2\alpha, 1\tau', ...
%                  'SF_C 2\alpha, 2\tau', ...
%                   ...
%                  'SF_F 1\alpha, 1\tau', ...
%                  'SF_F 1\alpha, 2\tau', ...
%                  'SF_F 2\alpha, 1\tau', ...
%                  'SF_F 2\alpha, 2\tau', ...
%                   ...
%                  'Trans 1\alpha, 1\tau', ...
%                  'Trans 1\alpha, 2\tau', ...
%                  'Trans 2\alpha, 1\tau', ...
%                  'Trans 2\alpha, 2\tau'};


sim_file_names = {'SF_compounds_2Alpha_2Tau_sim_data.mat', ...
                  ...
                  'SF_features_2Alpha_2Tau_sim_data.mat', ...
                   ...
                  'SF_features_transfer_choice_mixing_2Alpha_2Tau_sim_data.mat'};

model_names = {'SF_C 2\alpha, 2\tau', ...
                  ...
                 'SF_F 2\alpha, 2\tau', ...
                  ...
                 'Transfer 2\alpha, 2\tau'};

Labels_prior = {'4-cycle prior - 4-cycle',  '6-cycle prior - 4-line', '4-cycle prior - 6-line', '6-cycle prior - 6-cycle'};
Labels_trans = {'4-cycle prior - 4-cycle', '6-cycle prior - 4-cycle', '4-cycle prior - 6-cycle', '6-cycle prior - 6-cycle'};


for k = 1:size(model_names,2)

    if k == 1
       h = figure
       all_sims_plot = tiledlayout(3,2, 'TileSpacing', 'loose', 'Padding', 'loose');
    end


    %load data
    clear all_sim_data
    all_sim_data = load([res_dir, filesep, sim_file_names{k}]);
    all_sim_data = all_sim_data.sim_data;

    results = all_sim_data.results;
    res_to_plot = all_sim_data.res_to_plot;
    exp_res_to_plot = all_sim_data.exp_res_to_plot;
    inf_to_plot = all_sim_data.inf_to_plot;


    %transfer learning - experience probes
    nexttile
    avg_bar_scatter_simulations(results, exp_res_to_plot, Labels_trans, cols, [5 6], 2, ylim_lower, 4)
    ylabel('Probability Correct','FontSize',20,'FontName', 'Arial')

    if sum(k == [1 2]) > 0
        xticks([]);
        xticklabels({});
    end

    %transfer learning - inference probes
    nexttile
    avg_bar_scatter_simulations(results, inf_to_plot, Labels_trans, cols, [5 6], 2, ylim_lower, 4)

    ylabel('');
    if sum(k == [1 2]) > 0
        xticks([]);
        xticklabels({});
    end

    if k == 3
        screen_size = get(0, 'ScreenSize');
        origSize = get(h, 'Position'); % grab original on screen size
        set(h, 'Position', [0 0 screen_size(3)/1.75 screen_size(4) ] ); %set to scren size
        set(h,'PaperPositionMode','auto') %set paper pos for printing
        % saveas(h, 'fig_trajectories.png') % save figure
        ax = h;
        exportgraphics(ax,[plotpath 'fig_pcorrect_2a2t_models.png'],'Resolution',1200)
        set(h,'Position', origSize) %set back to original dimensions
    end
end






% lab_counter = 0;
% mod_counter = 0;

% for k = 1:size(model_names,2)
% 
%     if k == 1 || k == 5 || k == 9
%        h = figure
%        all_sims_plot = tiledlayout(2,8, 'TileSpacing', 'tight', 'Padding', 'tight');
%        lab_counter = 1;
%     end
% 
%     if sum(k == 2:4:size(model_names,2)) > 0
%        lab_counter = 1;
%     end
% 
%     if sum(k == [3 4 7 8 11 12]) > 0
%         mod_counter = 1;
%     end
% 
%     %load data
%     clear all_sim_data
%     all_sim_data = load([res_dir, filesep, sim_file_names{k}]);
%     all_sim_data = all_sim_data.sim_data;
% 
%     results = all_sim_data.results;
%     res_to_plot = all_sim_data.res_to_plot;
%     exp_res_to_plot = all_sim_data.exp_res_to_plot;
%     inf_to_plot = all_sim_data.inf_to_plot;
% 
%     nexttile
%     %prior learning - experience probes
%     avg_bar_scatter_simulations(results, exp_res_to_plot, Labels_prior, cols, [2 3], 1, ylim_lower, 4)
% 
%     if sum(k == 1:2:size(model_names,2)) > 0
%         ylabel('Probability Correct','FontSize',20,'FontName', 'Arial')
%     end
% 
%     if lab_counter == 1
%         title({'Prior: Exp', model_names{k}},'FontSize',20,'FontName', 'Arial');
%     end
% 
%     if mod_counter == 1
%        title(model_names{k},'FontSize',20,'FontName', 'Arial');
%     end
% 
%     if sum(k == [1 2 5 6 9 10]) > 0
%         xticks([]);
%         xticklabels({});
%     end
% 
% 
%     %prior learning - inference probes
%     nexttile
%     avg_bar_scatter_simulations(results, inf_to_plot, Labels_prior, cols, [2 3], 1, ylim_lower, 4)
%     if lab_counter == 1
%        title({'Prior: Inf', model_names{k}},'FontSize',20,'FontName', 'Arial');
%     end
% 
%     if mod_counter == 1
%        title(model_names{k},'FontSize',20,'FontName', 'Arial');
%     end
% 
%     if sum(k == [1 2 5 6 9 10]) > 0
%         xticks([]);
%         xticklabels({});
%     end
% 
%     %transfer learning - experience probes
%     nexttile
%     avg_bar_scatter_simulations(results, exp_res_to_plot, Labels_trans, cols, [5 6], 2, ylim_lower, 4)
%     if sum(k == 1:2:size(model_names,2)) > 0
%         ylabel('Probability Correct','FontSize',20,'FontName', 'Arial')
%     end
% 
%     if lab_counter == 1
%        title({'Transfer: Exp', model_names{k}},'FontSize',20,'FontName', 'Arial');
%     end
% 
% 
%     if mod_counter == 1
%        title(model_names{k},'FontSize',20,'FontName', 'Arial');
%     end
% 
%     ylabel('');
%     if sum(k == [1 2 5 6 9 10]) > 0
%         xticks([]);
%         xticklabels({});
%     end
% 
%     %transfer learning - inference probes
%     nexttile
%     avg_bar_scatter_simulations(results, inf_to_plot, Labels_trans, cols, [5 6], 2, ylim_lower, 4)
%     if lab_counter == 1
%        title({'Transfer: Inf', model_names{k}},'FontSize',20,'FontName', 'Arial');
%     end
% 
%     if mod_counter == 1
%        title(model_names{k},'FontSize',20,'FontName', 'Arial');
%     end
% 
%     ylabel('');
%     if sum(k == [1 2 5 6 9 10]) > 0
%         xticks([]);
%         xticklabels({});
%     end
% 
%     lab_counter = 0;
%     mod_counter = 0;
% 
%     % if sum(k == 4:4:12) > 0
%     %     screen_size = get(0, 'ScreenSize');
%     %     origSize = get(h, 'Position'); % grab original on screen size
%     %     set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size
%     %     set(h,'PaperPositionMode','auto') %set paper pos for printing
%     %     % saveas(h, 'fig_trajectories.png') % save figure
%     %     ax = h;
%     %     exportgraphics(ax,[plotpath 'fig_pcorrect_' model_names{k}(1:5) '.png'],'Resolution',300)
%     %     set(h,'Position', origSize) %set back to original dimensions
%     % end
% end







