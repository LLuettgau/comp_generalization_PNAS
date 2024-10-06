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
% ylabel('Probability Correct')
set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
ybounds = ylim;
yline(0.5,'--','linewidth',2.5)
%     set(gca,'YTick',[0 .25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'YTick',y_ticks, 'FontSize',26,'FontName', 'Arial');
set(gca,'TickDir','out')
set(gca,'xtick',1:4)
%set(gca,'XTickLabel', LabelsCS, 'FontSize',26,'FontName', 'Arial');
set(gca,'XTickLabelRotation', 45);
set(gcf,'color','w');
set(gca,'ycolor',cols.k)
set(gca,'xcolor',cols.k)

if mod_num == 4
    set(gca,'XTickLabel', LabelsCS, 'FontSize',26,'FontName', 'Arial');
    
    ticklabels_new = cell(size(LabelsCS));
    for i = 1:length(LabelsCS)
        ticklabels_new{i} = ['\color{black} ' LabelsCS{i}];
    end
    % set the tick labels
    set(gca, 'XTickLabel', ticklabels_new,'FontSize',20,'FontName', 'Arial');
end


% prepend a color for each tick label
LabelsY = get(gca,'YTickLabel');
ticklabels_ynew = cell(size(LabelsY));
for i = 1:length(LabelsY)
    ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
end
% set the tick labels
set(gca, 'YTickLabel', ticklabels_ynew,'FontSize',20,'FontName', 'Arial');  
%title('Overview','FontSize',25,'FontName', 'Arial');
