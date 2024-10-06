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


%% average (+data point) plots

%Overview
%plot p_correct 
positions = [1 2 3 4 5];   
pos = [positions-.275; ...
       positions+.275];
size_vec = [1 2 3 4 5];


markersize_overall = 50;
%prep CI
ci = 0.95;
alpha = 1-ci;

n = size(res_to_plot,1);
T_multiplier = tinv(1-alpha/2, n-1);

ci95 = T_multiplier*nanstd(res_to_plot(:,1:end))/sqrt(n);


g = figure('visible','on');
subplot(1,3,1)
b = bar(1:4,nanmean(res_to_plot(:,1:4)),'FaceColor','flat', 'BarWidth', 0.6, 'linewidth', 2.5, 'FaceAlpha',.5); hold all;
b.CData(1,:) = cols.dred;
b.CData(2,:) = cols.lyellow;
b.CData(3,:) = cols.medorange;
b.CData(4,:) = cols.cream;
errorbar(positions(1),nanmean(res_to_plot(:,1)), ci95(1), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
errorbar(positions(2),nanmean(res_to_plot(:,2)), ci95(2), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
errorbar(positions(3),nanmean(res_to_plot(:,3)), ci95(3), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
errorbar(positions(4),nanmean(res_to_plot(:,4)), ci95(4), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
%set(findobj(gca,'type','line'),'linew',2)
LabelsCS ={'Prior: Experience Probes', 'Prior: Inference Probes ', 'Transfer: Experience Probes', 'Transfer: Inference Probes'};
%ylim([0 1]);  
ylim([ylim_lower 1]);
xlim([0.5 4.5]);
box off
ylabel('Probability Correct')
set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
ybounds = ylim;
yline(0.5,'--','linewidth',2.5)
%     set(gca,'YTick',[0 .25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'YTick',y_ticks, 'FontSize',26,'FontName', 'Arial');
set(gca,'TickDir','out')
set(gca,'xtick',1:4)
set(gca,'XTickLabel', LabelsCS, 'FontSize',26,'FontName', 'Arial');
set(gca,'XTickLabelRotation', 45);
set(gcf,'color','w');
set(gca,'ycolor',cols.k)
set(gca,'xcolor',cols.k)
ticklabels_new = cell(size(LabelsCS));
for i = 1:length(LabelsCS)
    ticklabels_new{i} = ['\color{black} ' LabelsCS{i}];
end
% set the tick labels
set(gca, 'XTickLabel', ticklabels_new,'FontSize',20,'FontName', 'Arial');
% prepend a color for each tick label
LabelsY = get(gca,'YTickLabel');
ticklabels_ynew = cell(size(LabelsY));
for i = 1:length(LabelsY)
    ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
end
% set the tick labels
set(gca, 'YTickLabel', ticklabels_ynew,'FontSize',20,'FontName', 'Arial');  
%title('Overview','FontSize',25,'FontName', 'Arial');


%transfer learning - experience probes
Labels_trans = {'4-cycle prior - 4-cycle', '4-cycle prior - 6-cycle', '6-cycle prior - 4-cycle', '6-cycle prior - 6-cycle'};
subplot(1,3,2)
avg_bar_scatter_simulations(results, exp_res_to_plot, Labels_trans, cols, [5 6], 2, ylim_lower, 4)
%title('Transfer: Experience Probes','FontSize',22,'FontName', 'Arial');
%ylabel('Probability Correct','FontSize',20,'FontName', 'Arial')

%transfer learning - inference probes
subplot(1,3,3)
avg_bar_scatter_simulations(results, inf_to_plot, Labels_trans, cols, [5 6], 2, ylim_lower, 4)
%title('Transfer: Inference Probes','FontSize',22,'FontName', 'Arial');




%% broken up for conditions and trial types

Labels_prior = {'4-cycle prior - 4-cycle',  '6-cycle prior - 4-line', '4-cycle prior - 6-line', '6-cycle prior - 6-cycle'};
Labels_trans = {'4-cycle prior - 4-cycle', '6-cycle prior - 4-cycle', '4-cycle prior - 6-cycle', '6-cycle prior - 6-cycle'};

f = figure('visible','on');
%prior learning - experience probes
subplot(1,4,1)
avg_bar_scatter_simulations(results, exp_res_to_plot, Labels_prior, cols, [2 3], 1, ylim_lower, 4)
%title('Prior: Experience Probes','FontSize',22,'FontName', 'Arial');
ylabel('Probability Correct','FontSize',20,'FontName', 'Arial')


%prior learning - inference probes
subplot(1,4,2)
avg_bar_scatter_simulations(results, inf_to_plot, Labels_prior, cols, [2 3], 1, ylim_lower, 4)
%title('Prior: Inference Probes','FontSize',22,'FontName', 'Arial');


% h = figure('visible','on');
%transfer learning - experience probes
subplot(1,4,3)
avg_bar_scatter_simulations(results, exp_res_to_plot, Labels_trans, cols, [5 6], 2, ylim_lower, 4)
%title('Transfer: Experience Probes','FontSize',22,'FontName', 'Arial');
% ylabel('Probability Correct','FontSize',20,'FontName', 'Arial')

%transfer learning - inference probes
subplot(1,4,4)
avg_bar_scatter_simulations(results, inf_to_plot, Labels_trans, cols, [5 6], 2, ylim_lower, 4)
%title('Transfer: Inference Probes','FontSize',22,'FontName', 'Arial');

