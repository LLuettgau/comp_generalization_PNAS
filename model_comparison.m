%%model comparison using WAIC
clear all; close all; clc;


%% Define model labels for plots
addpath(genpath(['.', filesep, 'functions']))
modeldir = ['.', filesep, 'fitted_parameters/'];
datadir = ['.', filesep, 'data/'];
plotpath = ['.', filesep, 'fitted_parameters/model_comparison'];

models = {'fitted_parameters_MLE_SF_compounds_1Alpha_1Tau_rand_inival.mat', ...
          'fitted_parameters_MLE_SF_compounds_1Alpha_2Tau_rand_inival.mat', ...
          'fitted_parameters_MLE_SF_compounds_2Alpha_1Tau_rand_inival.mat', ...
          'fitted_parameters_MLE_SF_compounds_2Alpha_2Tau_rand_inival.mat', ...
          ...
          'fitted_parameters_MLE_SF_features_1Alpha_1Tau_rand_inival.mat', ...
          'fitted_parameters_MLE_SF_features_1Alpha_2Tau_rand_inival.mat', ...
          'fitted_parameters_MLE_SF_features_2Alpha_1Tau_rand_inival.mat', ...
          'fitted_parameters_MLE_SF_features_2Alpha_2Tau_rand_inival.mat', ...
          ...
          'fitted_parameters_MLE_SF_features_transfer_1Alpha_1Tau_1choiceMix_fitted_inival.mat', ...
          'fitted_parameters_MLE_SF_features_transfer_1Alpha_2Tau_1choiceMix_fitted_inival.mat', ...
          'fitted_parameters_MLE_SF_features_transfer_2Alpha_1Tau_1choiceMix_fitted_inival.mat', ...
          'fitted_parameters_MLE_SF_features_transfer_2Alpha_2Tau_1choiceMix_fitted_inival.mat'};

models_subset = {'fitted_parameters_MLE_SF_compounds_2Alpha_2Tau_rand_inival.mat', ...
          'fitted_parameters_MLE_SF_features_2Alpha_2Tau_rand_inival.mat', ...
          'fitted_parameters_MLE_SF_features_transfer_2Alpha_2Tau_1choiceMix_fitted_inival.mat'};

Labels_models_tmp ={'Compounds 1\alpha, 1\tau', ...
                 'Compounds 1\alpha, 2\tau', ...
                 'Compounds 2\alpha, 1\tau', ...
                 'Compounds 2\alpha, 2\tau', ...
                  ...
                 'Features 1\alpha, 1\tau', ...
                 'Features 1\alpha, 2\tau', ...
                 'Features 2\alpha, 1\tau', ...
                 'Features 2\alpha, 2\tau', ...
                  ...
                 'Transfer 1\alpha, 1\tau', ...
                 'Transfer 1\alpha, 2\tau', ...
                 'Transfer 2\alpha, 1\tau', ...
                 'Transfer 2\alpha, 2\tau'}; %...


all_data_tmp = load([datadir, filesep, 'all_data_agg.mat']);
exp_correct = all_data_tmp.data_all(:,1:7);
inf_correct = [all_data_tmp.data_all(:,1) all_data_tmp.data_all(:,9:end)];

%define font size for model labels in plots
font_size = 22;

%% get subIDs
subIDs_models_tmp = load([modeldir, 'subIDs_cell.mat']);
subIDs_models = subIDs_models_tmp.subject_ID_cell{1,1}(:);

conditions_models = subIDs_models_tmp.subject_ID_cell{1,2}(:);
scheds_models = subIDs_models_tmp.subject_ID_cell{1,3}(:);
scheds_models_unique = unique(scheds_models);

for i = 1:size(scheds_models,1)
    scheds_models_num(i,1) = find(strcmp(scheds_models(i),scheds_models_unique));
end

neglogliks = [];
all_params = struct();

for j = 1:size(models,2)
    data = importdata([modeldir, models{j}]);
    
    if j == 1
        subIDs_sorted = data.subIDs;
        ntr = data.nsamples;
    end

    field_name = genvarname(models{j}(1:end-4)); 
    all_params.(field_name) = data.params;

    neglogliks(:,j) = data.LLE;

    %calculate each model's WAIC
    for i = 1:size(data.params,1)
        % bic(i,j) = 2*(data.LLE(i))+size(data.params,2)*log(ntr(i)); 
       [total,se,pointwise] = waic(-1 * data.LLE(i));
        waic_all(i,j) = total.waic;
    end
end

%% alpha comparison
%compounds
mean(all_params.fitted_parameters_MLE_SF_compounds_1Alpha_1Tau_rand_inival(:,1))
mean(all_params.fitted_parameters_MLE_SF_compounds_2Alpha_1Tau_rand_inival(:,1:2))
mean(mean(all_params.fitted_parameters_MLE_SF_compounds_2Alpha_1Tau_rand_inival(:,1:2)))

mean(all_params.fitted_parameters_MLE_SF_compounds_1Alpha_2Tau_rand_inival(:,1))
mean(all_params.fitted_parameters_MLE_SF_compounds_2Alpha_2Tau_rand_inival(:,1:2))
mean(mean(all_params.fitted_parameters_MLE_SF_compounds_2Alpha_2Tau_rand_inival(:,1:2)))

%features
mean(all_params.fitted_parameters_MLE_SF_features_1Alpha_1Tau_rand_inival(:,1))
mean(all_params.fitted_parameters_MLE_SF_features_2Alpha_1Tau_rand_inival(:,1:2))
mean(mean(all_params.fitted_parameters_MLE_SF_features_2Alpha_1Tau_rand_inival(:,1:2)))

mean(all_params.fitted_parameters_MLE_SF_features_1Alpha_2Tau_rand_inival(:,1))
mean(all_params.fitted_parameters_MLE_SF_features_2Alpha_2Tau_rand_inival(:,1:2))
mean(mean(all_params.fitted_parameters_MLE_SF_features_2Alpha_2Tau_rand_inival(:,1:2)))

%features transfer
mean(all_params.fitted_parameters_MLE_SF_features_transfer_1Alpha_1Tau_1choiceM(:,1))
mean(all_params.fitted_parameters_MLE_SF_features_transfer_2Alpha_1Tau_1choiceM(:,1:2))
mean(mean(all_params.fitted_parameters_MLE_SF_features_transfer_2Alpha_1Tau_1choiceM(:,1:2)))

mean(all_params.fitted_parameters_MLE_SF_features_transfer_1Alpha_2Tau_1choiceM(:,1))
mean(all_params.fitted_parameters_MLE_SF_features_transfer_2Alpha_2Tau_1choiceM(:,1:2))
mean(mean(all_params.fitted_parameters_MLE_SF_features_transfer_2Alpha_2Tau_1choiceM(:,1:2)))


%% plot parameter estimate distributions across models

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

alpha_color = cols.dred;
tau_color = cols.b;
omega_color = cols.dgrey;

bin_num_alpha = 50;
bin_num_tau = 241;
bin_num_omega = 50;

xlims = [-0.01 1; -1 100];
ylims = [0 150; 0 160; 0 15];

h = figure
all_params_plot = tiledlayout(3,5, 'TileSpacing', 'tight', 'Padding', 'tight');
set(gcf, 'Color', 'w'); % 'w' stands for white
ylabel(all_params_plot,'#Subjects', 'FontSize',26,'FontName', 'Arial')

for i = 1:size(models_subset,2)
    
    field_name = genvarname(models_subset{i}(1:end-4));

    %Compounds / Features
    if i == 1 || i == 2


        nexttile    
        histogram(all_params.(field_name)(:,1), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
        ylim(ylims(1,:));
        xlim(xlims(1,:));
        set(gca,'YTick',ylims(1,1):50:ylims(1,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'YTickLabel', ylims(1,1):50:ylims(1,2),'FontSize',18,'FontName', 'Arial');
        box off 
        set(gca, 'TickDir', 'out');
        xticks([]);
        xticklabels({});
        if i == 1
            title('\alpha_1', 'FontSize',22,'FontName', 'Arial');
        end

        nexttile    
        histogram(all_params.(field_name)(:,2), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
        ylim(ylims(1,:));
        xlim(xlims(1,:));
        set(gca,'YTick',ylims(1,1):50:ylims(1,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'YTickLabel', ylims(1,1):50:ylims(1,2),'FontSize',18,'FontName', 'Arial');
        box off 
        set(gca, 'TickDir', 'out');
        xticks([]);
        xticklabels({});
        if i == 1
            title('\alpha_2', 'FontSize',22,'FontName', 'Arial');
        end

        nexttile    
        histogram(all_params.(field_name)(:,3), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
        ylim(ylims(2,:));
        xlim(xlims(2,:));
        box off 
        set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
        set(gca, 'TickDir', 'out');
        xticks([]);
        xticklabels({});
        if i == 1
            title('\tau_1', 'FontSize',22,'FontName', 'Arial');
        end

        nexttile    
        histogram(all_params.(field_name)(:,4), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
        ylim(ylims(2,:));
        xlim(xlims(2,:));
        box off 
        set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
        set(gca, 'TickDir', 'out');
        xticks([]);
        xticklabels({});
        if i == 1
            title('\tau_2', 'FontSize',22,'FontName', 'Arial');
        end

        nexttile
        if i == 1
            title('\omega', 'FontSize',22,'FontName', 'Arial');
        end
        axis off;

    else

    %Transfer

        nexttile    
        histogram(all_params.(field_name)(:,1), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
        ylim(ylims(1,:));
        xlim(xlims(1,:));
        box off 
        set(gca,'XTick',0:.25:xlims(1,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'XTickLabel', 0:.25:xlims(1,2),'FontSize',18,'FontName', 'Arial');
        set(gca, 'TickDir', 'out');
        xtickangle(45);

        nexttile    
        histogram(all_params.(field_name)(:,2), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
        ylim(ylims(1,:));
        xlim(xlims(1,:));
        box off 
        set(gca,'XTick',0:.25:xlims(1,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'XTickLabel', 0:.25:xlims(1,2),'FontSize',18,'FontName', 'Arial');
        set(gca, 'TickDir', 'out');
        xtickangle(45);

        nexttile    
        histogram(all_params.(field_name)(:,3), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
        ylim(ylims(2,:));
        xlim(xlims(2,:));
        box off 
        set(gca,'XTick',0:25:xlims(2,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'XTickLabel', 0:25:xlims(2,2),'FontSize',18,'FontName', 'Arial');
        set(gca,'YTick',ylims(2,1):60:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'YTickLabel', ylims(2,1):60:ylims(2,2),'FontSize',18,'FontName', 'Arial');
        set(gca, 'TickDir', 'out');
        xtickangle(45);

        nexttile    
        histogram(all_params.(field_name)(:,4), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
        ylim(ylims(2,:));
        xlim(xlims(2,:));
        box off 
        set(gca,'XTick',0:25:xlims(2,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'XTickLabel', 0:25:xlims(2,2),'FontSize',18,'FontName', 'Arial');
        set(gca,'YTick',ylims(2,1):60:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'YTickLabel', ylims(2,1):60:ylims(2,2),'FontSize',18,'FontName', 'Arial');
        set(gca, 'TickDir', 'out');
        xtickangle(45);

        nexttile    
        histogram(all_params.(field_name)(:,5), 'FaceColor', omega_color, 'EdgeColor', omega_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_omega);
        ylim(ylims(3,:));
        xlim(xlims(1,:));
        box off 
        set(gca,'XTick',0:.25:xlims(1,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'XTickLabel',0:.25:xlims(1,2),'FontSize',18,'FontName', 'Arial');
        set(gca,'YTick',ylims(3,1):5:ylims(3,2), 'FontSize',18,'FontName', 'Arial');
        set(gca, 'YTickLabel', ylims(3,1):5:ylims(3,2),'FontSize',18,'FontName', 'Arial');
        set(gca, 'TickDir', 'out');
        xtickangle(45);

    end


end


% screen_size = get(0, 'ScreenSize');
% origSize = get(h, 'Position'); % grab original on screen size
% set(h, 'Position', [0 0 screen_size(3)/1.25 screen_size(4)/1.5 ] ); %set to scren size
% set(h,'PaperPositionMode','auto') %set paper pos for printing
% ax = h;
% exportgraphics(ax,[plotpath, filesep, 'param_estimate_dists_best_versions.png'],'Resolution',1200)
% set(h,'Position', origSize) %set back to original dimensions
% 
% 
% %plot all parameter distributions
% for i = 1:size(models,2)
% 
%     field_name = genvarname(models{i}(1:end-4));
% 
%     if i == 1 || i == 5      
% 
%         num_cols = 4;
% 
%         h = figure
%         all_params_plot = tiledlayout(4,num_cols, 'TileSpacing', 'loose', 'Padding', 'tight');
%         set(gcf, 'Color', 'w'); % 'w' stands for white
%         ylabel(all_params_plot,'#Subjects', 'FontSize',26,'FontName', 'Arial')
%     elseif i == 9
% 
%         num_cols = 5;
% 
%         h = figure
%         all_params_plot = tiledlayout(4,num_cols, 'TileSpacing', 'loose', 'Padding', 'loose');
%         set(gcf, 'Color', 'w'); % 'w' stands for white
%         ylabel(all_params_plot,'#Subjects', 'FontSize',26,'FontName', 'Arial')
% 
%     end
% 
%     %1alpha/1tau
%     if size(all_params.(field_name),2) == 2
%         nexttile;
%         histogram(all_params.(field_name)(:,1), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         if i == 1 || i == 5 
%             title('\alpha_1', 'FontSize',28,'FontName', 'Arial');
%         end
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(1,1):50:ylims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(1,1):50:ylims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile
%         title('\alpha_2', 'FontSize',22,'FontName', 'Arial');
%         axis off;
% 
%         nexttile
%         histogram(all_params.(field_name)(:,2), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         if i == 1 || i == 5 
%             title('\tau_1', 'FontSize',28,'FontName', 'Arial');
%         end
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile
%         title('\tau_2', 'FontSize',22,'FontName', 'Arial');
%         axis off;
% 
%     %1alpha/2tau
%     elseif size(all_params.(field_name),2) == 3 && sum(i == [2 6]) > 0
%         nexttile    
%         histogram(all_params.(field_name)(:,1), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(1,1):50:ylims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(1,1):50:ylims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile
%         axis off;
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,2), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,3), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
% 
%     %2alpha/1tau
%     elseif size(all_params.(field_name),2) == 3 && sum(i == [3 7]) > 0
%         nexttile    
%         histogram(all_params.(field_name)(:,1), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(1,1):50:ylims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(1,1):50:ylims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,2), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(2,1):50:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):50:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,3), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         axis off;
% 
%     %2alpha/2tau
%     elseif size(all_params.(field_name),2) == 4 && sum(i == [10 11]) == 0
%         nexttile    
%         histogram(all_params.(field_name)(:,1), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         set(gca,'XTick',0:.25:xlims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'XTickLabel', 0:.25:xlims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
%         xtickangle(45);
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,2), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         set(gca,'XTick',0:.25:xlims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'XTickLabel', 0:.25:xlims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
%         xtickangle(45);
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,3), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         set(gca,'XTick',0:25:xlims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'XTickLabel', 0:25:xlims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
%         xtickangle(45);
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,4), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         set(gca,'XTick',0:25:xlims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'XTickLabel', 0:25:xlims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
%         xtickangle(45);
% 
%     %1alpha/1tau/1mix
%     elseif sum(i == 9) > 0
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,1), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         title('\alpha_1', 'FontSize',28,'FontName', 'Arial');
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(1,1):50:ylims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(1,1):50:ylims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile
%         title('\alpha_2', 'FontSize',22,'FontName', 'Arial');
%         axis off;
% 
%         nexttile
%         histogram(all_params.(field_name)(:,2), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         title('\tau_1', 'FontSize',28,'FontName', 'Arial');
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile
%         title('\tau_2', 'FontSize',22,'FontName', 'Arial');
%         axis off;
% 
%         nexttile
%         histogram(all_params.(field_name)(:,3), 'FaceColor', omega_color, 'EdgeColor', omega_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_omega);
%         ylim(ylims(3,:));
%         xlim(xlims(1,:));
%         box off 
%         title('\omega', 'FontSize',28,'FontName', 'Arial');
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(3,1):5:ylims(3,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(3,1):5:ylims(3,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%     %1alpha/2tau/1mix
%     elseif sum(i == 10) > 0
%         nexttile    
%         histogram(all_params.(field_name)(:,1), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(1,1):50:ylims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(1,1):50:ylims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         axis off;
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,2), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,3), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,4), 'FaceColor', omega_color, 'EdgeColor', omega_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_omega);
%         ylim(ylims(3,:));
%         xlim(xlims(1,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(3,1):5:ylims(3,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(3,1):5:ylims(3,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%     %2alpha/1tau/1mix
%     elseif sum(i == 11) > 0
%         nexttile    
%         histogram(all_params.(field_name)(:,1), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(1,1):50:ylims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(1,1):50:ylims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,2), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(1,1):50:ylims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(1,1):50:ylims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,3), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(2,1):40:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):40:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         axis off;
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,4), 'FaceColor', omega_color, 'EdgeColor', omega_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_omega);
%         ylim(ylims(3,:));
%         xlim(xlims(1,:));
%         box off 
%         xticks([]);
%         xticklabels({});
%         set(gca,'YTick',ylims(3,1):5:ylims(3,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(3,1):5:ylims(3,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%     %2alpha/2tau/1mix
%     elseif sum(i == 12) > 0
%         nexttile    
%         histogram(all_params.(field_name)(:,1), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         set(gca,'XTick',0:.25:xlims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'XTickLabel', 0:.25:xlims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,2), 'FaceColor', alpha_color, 'EdgeColor', alpha_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_alpha);
%         ylim(ylims(1,:));
%         xlim(xlims(1,:));
%         box off 
%         set(gca,'XTick',0:.25:xlims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'XTickLabel', 0:.25:xlims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,3), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         set(gca,'XTick',0:25:xlims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'XTickLabel', 0:25:xlims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca,'YTick',ylims(2,1):60:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):60:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,4), 'FaceColor', tau_color, 'EdgeColor', tau_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_tau);
%         ylim(ylims(2,:));
%         xlim(xlims(2,:));
%         box off 
%         set(gca,'XTick',0:25:xlims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'XTickLabel', 0:25:xlims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca,'YTick',ylims(2,1):60:ylims(2,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(2,1):60:ylims(2,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%         nexttile    
%         histogram(all_params.(field_name)(:,5), 'FaceColor', omega_color, 'EdgeColor', omega_color, 'LineWidth', 1, 'FaceAlpha', 0.75, 'NumBins', bin_num_omega);
%         ylim(ylims(3,:));
%         xlim(xlims(1,:));
%         box off 
%         set(gca,'XTick',0:.25:xlims(1,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'XTickLabel',0:.25:xlims(1,2),'FontSize',18,'FontName', 'Arial');
%         set(gca,'YTick',ylims(3,1):5:ylims(3,2), 'FontSize',18,'FontName', 'Arial');
%         set(gca, 'YTickLabel', ylims(3,1):5:ylims(3,2),'FontSize',18,'FontName', 'Arial');
%         set(gca, 'TickDir', 'out');
% 
%     end
% 
%     %  if sum(i == [4 8 12]) > 0
%     %     screen_size = get(0, 'ScreenSize');
%     %     origSize = get(h, 'Position'); % grab original on screen size
%     %     set(h, 'Position', [0 0 screen_size(3)/1.5 screen_size(4)/2 ] ); %set to scren size
%     %     set(h,'PaperPositionMode','auto') %set paper pos for printing
%     %     ax = h;
%     %     exportgraphics(ax,[plotpath, filesep, 'param_estimate_dists_' field_name(23:39) '.png'],'Resolution',1200)
%     %     set(h,'Position', origSize) %set back to original dimensions
%     % end
% 
% end



%% put WAIC and LLE to tables
waic_table = array2table(waic_all, 'VariableNames', Labels_models_tmp(1:end));
mean_waic_table = mean(waic_table);

LLE_table = array2table(neglogliks, 'VariableNames', Labels_models_tmp(1:end));

nanmedian(data.params(conditions_models == 1,1))
nanmedian(data.params(conditions_models == 2,1))
nanmean(data.params(conditions_models == 1,1))
nanmean(data.params(conditions_models == 2,1))

nanmedian(data.params(conditions_models == 1,2))
nanmedian(data.params(conditions_models == 2,2))

nanmedian(data.params(conditions_models == 1,3))
nanmedian(data.params(conditions_models == 2,3))

nanmedian(data.params(conditions_models == 1,4))
nanmedian(data.params(conditions_models == 2,4))


figure; 
subplot(1,2,1)
hist(data.params(conditions_models == 1,1),50)
subplot(1,2,2)
hist(data.params(conditions_models == 2,1),50)

figure; 
subplot(1,2,1)
hist(data.params(conditions_models == 1,2),50)
subplot(1,2,2)
hist(data.params(conditions_models == 2,2),50)

figure; 
subplot(1,2,1)
hist(data.params(conditions_models == 1,3),50)
subplot(1,2,2)
hist(data.params(conditions_models == 2,3),50)

figure; 
subplot(1,2,1)
hist(data.params(conditions_models == 1,4),50)
subplot(1,2,2)
hist(data.params(conditions_models == 2,4),50)


%correlation between fitted params
for i = 1:3
    for j = 1:3
        [r_params(i,j), p_params(i,j)] = corr(data.params(:,i), data.params(:,j), 'type', 'Pearson', 'rows', 'complete');
    end
end

xlabels = {'Alpha', 'Tau', 'Mixing'};

figure
h = heatmap(r_params); colorbar
h.XDisplayLabels = xlabels;
h.YDisplayLabels = xlabels;

%parameter and performance correlations
%correlation between fitted params
%mixing parameter
for i = 2:size(exp_correct,2)
    [r_corr_exp(i), p_corr_exp(i)] = corr(exp_correct(:,i), data.params(:,3), 'type', 'Pearson', 'rows', 'complete');
end

for i = 2:size(exp_correct,2)
    [r_corr_inf(i), p_corr_inf(i)] = corr(inf_correct(:,i), data.params(:,3), 'type', 'Pearson', 'rows', 'complete');
end


%% 
%compute WAIC for structural inference model
si_model_tmp = readtable([modeldir, 'fits_all_llmax.csv']);
si_model(:,1) = si_model_tmp.subj;
si_model(:,2) = si_model_tmp.subID;
si_model(:,3) = -1*si_model_tmp.llmax;
si_model(:,4) = si_model_tmp.bic;


for k = 1:size(si_model,1)
    if isinf(si_model(k,3))
        si_model(k,3) = 116.63;
        si_model(k,4) = 2*(116.63)+log(164);
        [total,se,pointwise] = waic(-1 * 116.63);
        si_model(k,5) = total.waic;
    else
       [total,se,pointwise] = waic(-1 * si_model(k,3));
       si_model(k,5) = total.waic;
    end
end


% waic_all = [waic_all si_model(:,5)];
% neglogliks_all = [neglogliks si_model(:,3)];

%sort labels and BIC scores from best to worst model 
waic_res_all = waic_all;
waic_tmp = waic_res_all;
[~,idx] = sort(mean(waic_tmp));

for k = 1:size(waic_tmp,2)
    waic_res_all(:,k) = waic_tmp(:,idx(k));
    Labels_waic(k) = Labels_models_tmp(idx(k));
end


%get model frequencies --> find column of subjects for who this WAIC
%for this model is lowest
modelfreq(1,:) = sum(waic_res_all(:,1:end) == min(waic_res_all(:,1:end),[],2));
rel_modelfreq = modelfreq / size(waic_res_all,1);

disp('mean WAIC')
mean(waic_res_all)
winnerWAIC = Labels_waic(mean(waic_res_all)==min(mean(waic_res_all)));
min(mean(waic_res_all))

%total of subjects for which transfer models are best
sum(modelfreq([1 2 4 6]))
sum(rel_modelfreq([1 2 4 6]))

%% Plot WAIC and model frequencies
% %prepare confidence interval
ci = 0.95 ;
alpha = 1-ci;

n = size(waic_res_all,1);
T_multiplier = tinv(1-alpha/2, n-1);

ci95 = T_multiplier*nanstd(waic_res_all)/sqrt(n);


% positions = 1:size(models,2)+1;   
positions = 1:size(models,2);   
pos = [positions-.2; ...
       positions+.2];
   
%% plot WAIC scores

h = figure
mod_comp_res = tiledlayout(3,1, 'TileSpacing', 'loose', 'Padding', 'tight');
nexttile

bar(1:size(waic_res_all,2), mean(waic_res_all), 'w', 'BarWidth', 0.6, 'linewidth', 2.5); hold all;
for k = 1:size(positions,2)
    errorbar(positions(k), mean(waic_res_all(:,k)), ci95(k), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
end
ylabel('WAIC', 'FontSize',30,'FontName', 'Arial')
ylim([170 210]);     
xtickangle(45)
hold off
box off
set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
set(findobj(gca,'type','line'),'linew',2)
set(gcf,'color','w');
set(gca,'TickDir','out')
set(gca,'xtick',1:size(waic_res_all,2))
set(gca,'XTickLabel', Labels_waic, 'FontSize',font_size,'FontName', 'Arial');
set(gca,'XTickLabelRotation', 45);
set(gcf,'color','w');
ticklabels_new = cell(size(Labels_waic));
for i = 1:length(Labels_waic)
    ticklabels_new{i} = ['\color{black} ' Labels_waic{i}];
end
set(gca,'YTick',170:10:210, 'FontSize',26,'FontName', 'Arial');
set(gca, 'YTickLabel', 170:10:210,'FontSize',20,'FontName', 'Arial');
xticks([]);
xticklabels({});

% screen_size = get(0, 'ScreenSize');
% origSize = get(h, 'Position'); % grab original on screen size
% set(h, 'Position', [0 0 screen_size(3)/1.5 screen_size(4)/1.5 ] ); %set to scren size
% set(h,'PaperPositionMode','auto') %set paper pos for printing
% ax = h;
% exportgraphics(ax,[plotpath, filesep, 'WAIC_scores.png'],'Resolution',300)
% set(h,'Position', origSize) %set back to original dimensions
% 


%% plot WAIC differences relative to best-fitting model

waic_res_all_diff = waic_res_all - waic_res_all(:,1);

% %prepare confidence interval
n = size(waic_res_all_diff,1);
T_multiplier = tinv(1-alpha/2, n-1);

ci95_diff = T_multiplier*nanstd(waic_res_all_diff)/sqrt(n);


% h = figure
nexttile

bar(1:size(waic_res_all_diff,2), mean(waic_res_all_diff), 'w', 'BarWidth', 0.6, 'linewidth', 2.5); hold all;
for k = 1:size(positions,2)
    errorbar(positions(k), mean(waic_res_all_diff(:,k)), ci95_diff(k), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
end
ylabel('\Delta WAIC', 'FontSize',30,'FontName', 'Arial')
ylim([0 30]);     
xtickangle(45)
hold off
box off
set(gcf,'color','w');
set(gca,'TickDir','out')
set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
set(findobj(gca,'type','line'),'linew',2)

set(gca,'xtick',1:size(waic_res_all_diff,2))
set(gca,'XTickLabel', Labels_waic, 'FontSize',font_size,'FontName', 'Arial');
set(gca,'XTickLabelRotation', 45);
set(gcf,'color','w');
ticklabels_new = cell(size(Labels_waic));
for i = 1:length(Labels_waic)
    ticklabels_new{i} = ['\color{black} ' Labels_waic{i}];
end
set(gca,'YTick',0:10:30, 'FontSize',26,'FontName', 'Arial');
set(gca, 'YTickLabel', 0:10:30,'FontSize',20,'FontName', 'Arial');
xticks([]);
xticklabels({});

% screen_size = get(0, 'ScreenSize');
% origSize = get(h, 'Position'); % grab original on screen size
% set(h, 'Position', [0 0 screen_size(3)/1.5 screen_size(4)/1.5 ] ); %set to scren size
% set(h,'PaperPositionMode','auto') %set paper pos for printing
% ax = h;
% exportgraphics(ax,[plotpath, filesep, 'WAIC_model_differences.png'],'Resolution',300)
% set(h,'Position', origSize) %set back to original dimensions


%% plot model frequencies
freq_data = modelfreq / sum(modelfreq);
freq_data_models(1) = sum(rel_modelfreq([1 2 4 6]));
freq_data_models(2) = sum(rel_modelfreq([3 5 7 8]));
freq_data_models(3) = sum(rel_modelfreq(9:12));
Labels_mod = {'Transfer', ...
              'Features', ...
              'Compounds'};

Labels = Labels_waic;
ylbl = 'WAIC Model Frequencies';

% h = figure
nexttile

bar(1:size(freq_data,2), freq_data, 'k', 'BarWidth', 0.6, 'linewidth', 2.5); hold all;
ylabel(ylbl, 'FontSize',25,'FontName', 'Arial')
ylim([0 0.60]);
xtickangle(45)
hold off
box off
set(gcf,'color','w');
set(gca,'TickDir','out')
set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
set(findobj(gca,'type','line'),'linew',2)
set(gca,'xtick',1:length(freq_data))
set(gca,'XTickLabel', Labels, 'FontSize',font_size,'FontName', 'Arial');
set(gca,'XTickLabelRotation', 45);
set(gcf,'color','w');
ticklabels_new = cell(size(Labels));
for i = 1:length(Labels)
    ticklabels_new{i} = ['\color{black} ' Labels{i}];
end
set(gca,'YTick',0:.2:0.6, 'FontSize',26,'FontName', 'Arial');
set(gca, 'YTickLabel', 0:.2:0.6,'FontSize',20,'FontName', 'Arial');
% Plotting the histogram
inset_ax = axes('Position',[0.4, .282, 1.01, 0.1]); 
box(inset_ax, 'on');
bar(inset_ax, 1:size(freq_data_models,2), freq_data_models, 'k', 'BarWidth', 0.6, 'linewidth', 2.5); hold all;
box off
set(gca,'YTick',0:.5:1, 'FontSize',26,'FontName', 'Arial');
set(gca, 'YTickLabel', 0:.5:1,'FontSize',20,'FontName', 'Arial');
set(gcf,'color','w');
set(gca,'TickDir','out')
set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
set(findobj(gca,'type','line'),'linew',2)
pbaspect([1 1.1 1])
xlim([-.01 3.5]);
set(gca,'XTickLabel', Labels_mod, 'FontSize',font_size-2,'FontName', 'Arial');
set(gca,'XTickLabelRotation', 45);
ticklabels_new = cell(size(Labels_mod));
for i = 1:length(Labels_mod)
    ticklabels_new{i} = ['\color{black} ' Labels_mod{i}];
end
% screen_size = get(0, 'ScreenSize');
% origSize = get(h, 'Position'); % grab original on screen size
% set(h, 'Position', [0 0 screen_size(3)/1.5 screen_size(4)/1.5 ] ); %set to scren size
% set(h,'PaperPositionMode','auto') %set paper pos for printing
% ax = h;
% exportgraphics(ax,[plotpath, filesep, 'WAIC_model_frequencies.png'],'Resolution',300)
% set(h,'Position', origSize) %set back to original dimensions

screen_size = get(0, 'ScreenSize');
origSize = get(h, 'Position'); % grab original on screen size
set(h, 'Position', [0 0 screen_size(3)/3.5 screen_size(4) ] ); %set to scren size
set(h,'PaperPositionMode','auto') %set paper pos for printing
ax = h;
exportgraphics(ax,[plotpath, filesep, 'WAIC_model_comparison_vert.png'],'Resolution',1200)
set(h,'Position', origSize) %set back to original dimensions

