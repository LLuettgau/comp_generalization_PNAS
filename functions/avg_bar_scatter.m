%Function to plot bar plots and scatter for condition averages

function [] = avg_bar_scatter(results, sum_data, markersize, Labels, cols, col_idx, phase, probe_type)

    %plot p_correct broken up for conditions
    if sum(results(:,8) == 1,1) > sum(results(:,8) == 2,1)
        data_to_plot = NaN(sum(results(:,8) == 1,1),4);
    else
        data_to_plot = NaN(sum(results(:,8) == 2,1),4);
    end
    %ordered by condition
%     data_to_plot(1:sum(results(:,8) == 1),1) = sum_data(results(:,8) == 1,col_idx(1));
%     data_to_plot(1:sum(results(:,8) == 1),2) = sum_data(results(:,8) == 1,col_idx(2));
%     data_to_plot(1:sum(results(:,8) == 2),3) = sum_data(results(:,8) == 2,col_idx(1));
%     data_to_plot(1:sum(results(:,8) == 2),4) = sum_data(results(:,8) == 2,col_idx(2));
    
    %ordered by trial type
    data_to_plot(1:sum(results(:,8) == 1),1) = sum_data(results(:,8) == 1,col_idx(1));
    data_to_plot(1:sum(results(:,8) == 2),2) = sum_data(results(:,8) == 2,col_idx(1));
    data_to_plot(1:sum(results(:,8) == 1),3) = sum_data(results(:,8) == 1,col_idx(2));
    data_to_plot(1:sum(results(:,8) == 2),4) = sum_data(results(:,8) == 2,col_idx(2));
    
    

    positions = [1 2 3 4 5];   
    pos = [positions-.3; ...
           positions+.3];
    size_vec = [1 2 3 4 5];

    % %prepare confidence interval
    ci = 0.95 ;
    alpha = 1-ci;

    n = size(data_to_plot,1);
    T_multiplier = tinv(1-alpha/2, n-1);

    ci95 = T_multiplier*nanstd(data_to_plot)/sqrt(n);

    bar(1:4,nanmean(data_to_plot), 'w', 'BarWidth', 0.6, 'linewidth', 2.5); hold all;
%     coords1 = linspace(pos(1,1), pos(2,1),length(data_to_plot(:,1)));
%     coords2 = linspace(pos(1,2), pos(2,2), length(data_to_plot(:,2)));
%     coords3 = linspace(pos(1,3), pos(2,3), length(data_to_plot(:,3)));
%     coords4 = linspace(pos(1,4), pos(2,4), length(data_to_plot(:,4)));

    %ordered by condition

% %     errorbar(positions(1),nanmean(data_to_plot(:,1)), ci95(1), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
% %     errorbar(positions(2),nanmean(data_to_plot(:,2)), ci95(2), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
% %     errorbar(positions(3),nanmean(data_to_plot(:,3)), ci95(3), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
% %     errorbar(positions(4),nanmean(data_to_plot(:,4)), ci95(4), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);     
%     set(findobj(gca,'type','line'),'linew',2)
%     
%     for i = 1:size(data_to_plot,1)
%         plot([coords1(i) coords2(i)], [data_to_plot(i,1) data_to_plot(i,2)], '-', 'LineWidth', 1, 'Color', cols.grey)
%         plot([coords3(i) coords4(i)], [data_to_plot(i,3) data_to_plot(i,4)], '-', 'LineWidth', 1, 'Color', cols.grey)
%     end
%     scatter(linspace(pos(1,1), pos(2,1),length(data_to_plot(:,1))), data_to_plot(:,1), markersize, 'o', 'MarkerFaceColor',  cols.y, 'MarkerEdgeColor', cols.y,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
%     scatter(linspace(pos(1,2), pos(2,2), length(data_to_plot(:,2))), data_to_plot(:,2), markersize, 'o', 'MarkerFaceColor',  cols.ytrend, 'MarkerEdgeColor', cols.ytrend,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
%     scatter(linspace(pos(1,3), pos(2,3), length(data_to_plot(:,3))), data_to_plot(:,3), markersize, 'o', 'MarkerFaceColor',  cols.b, 'MarkerEdgeColor', cols.b,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
%     scatter(linspace(pos(1,4), pos(2,4), length(data_to_plot(:,4))), data_to_plot(:,4), markersize, 'o', 'MarkerFaceColor',  cols.btrend, 'MarkerEdgeColor', cols.btrend,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
    
% if phase == 1
%     scatter(linspace(pos(1,1), pos(2,1),length(data_to_plot(:,1))), data_to_plot(:,1), markersize, 'o', 'MarkerFaceColor',  cols.y, 'MarkerEdgeColor', cols.y,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
%     scatter(linspace(pos(1,2), pos(2,2), length(data_to_plot(:,2))), data_to_plot(:,2), markersize, 'o', 'MarkerFaceColor',  cols.ytrend, 'MarkerEdgeColor', cols.ytrend,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
%     scatter(linspace(pos(1,3), pos(2,3), length(data_to_plot(:,3))), data_to_plot(:,3), markersize, 'o', 'MarkerFaceColor',  cols.b, 'MarkerEdgeColor', cols.b,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
%     scatter(linspace(pos(1,4), pos(2,4), length(data_to_plot(:,4))), data_to_plot(:,4), markersize, 'o', 'MarkerFaceColor',  cols.btrend, 'MarkerEdgeColor', cols.btrend,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
% else
%     scatter(linspace(pos(1,1), pos(2,1),length(data_to_plot(:,1))), data_to_plot(:,1), markersize, 'o', 'MarkerFaceColor',  cols.green, 'MarkerEdgeColor', cols.green,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
%     scatter(linspace(pos(1,2), pos(2,2), length(data_to_plot(:,2))), data_to_plot(:,2), markersize, 'o', 'MarkerFaceColor',  cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
%     scatter(linspace(pos(1,3), pos(2,3), length(data_to_plot(:,3))), data_to_plot(:,3), markersize, 'o', 'MarkerFaceColor',  cols.green, 'MarkerEdgeColor', cols.green,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
%     scatter(linspace(pos(1,4), pos(2,4), length(data_to_plot(:,4))), data_to_plot(:,4), markersize, 'o', 'MarkerFaceColor',  cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
% end

%ordered by trial type
if phase == 1
    scatter(linspace(pos(1,1), pos(2,1),length(data_to_plot(:,1))), data_to_plot(:,1), markersize, 'o', 'MarkerFaceColor',  cols.y, 'MarkerEdgeColor', cols.y,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
    scatter(linspace(pos(1,2), pos(2,2), length(data_to_plot(:,2))), data_to_plot(:,2), markersize, 'o', 'MarkerFaceColor',  cols.b, 'MarkerEdgeColor', cols.b,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
    scatter(linspace(pos(1,3), pos(2,3), length(data_to_plot(:,3))), data_to_plot(:,3), markersize, 'o', 'MarkerFaceColor',  cols.ytrend, 'MarkerEdgeColor', cols.ytrend,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
    scatter(linspace(pos(1,4), pos(2,4), length(data_to_plot(:,4))), data_to_plot(:,4), markersize, 'o', 'MarkerFaceColor',  cols.btrend, 'MarkerEdgeColor', cols.btrend,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
else
    scatter(linspace(pos(1,1), pos(2,1),length(data_to_plot(:,1))), data_to_plot(:,1), markersize, 'o', 'MarkerFaceColor',  cols.green, 'MarkerEdgeColor', cols.green,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
    scatter(linspace(pos(1,2), pos(2,2), length(data_to_plot(:,2))), data_to_plot(:,2), markersize, 'o', 'MarkerFaceColor',  cols.green, 'MarkerEdgeColor', cols.green,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
    scatter(linspace(pos(1,3), pos(2,3), length(data_to_plot(:,3))), data_to_plot(:,3), markersize, 'o', 'MarkerFaceColor',  cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
    scatter(linspace(pos(1,4), pos(2,4), length(data_to_plot(:,4))), data_to_plot(:,4), markersize, 'o', 'MarkerFaceColor',  cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
end




    errorbar(positions(1),nanmean(data_to_plot(:,1)), ci95(1), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
    errorbar(positions(2),nanmean(data_to_plot(:,2)), ci95(2), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
    errorbar(positions(3),nanmean(data_to_plot(:,3)), ci95(3), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
    errorbar(positions(4),nanmean(data_to_plot(:,4)), ci95(4), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);     
    
    ylim([0 1]);     
    %xlim([0 5]);
    xlim([0.5 4.5]);
    box off
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)        
    set(gca,'YTick',[0 .25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');

    ybounds = ylim;
    yline(0.5,'--') 
    %set(gca,'YTick',[0 .25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'xtick',1:4)
    set(gca,'XTickLabel', Labels, 'FontSize',20,'FontName', 'Arial');
    set(gca,'XTickLabelRotation', 45);
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
    
    if probe_type == 2
        %         set(gca,'YTick',[])
        %         ax1 = gca;
        %         ax1.YAxis.Visible = 'off';   % remove y-axis
        
        set(gca,'XTick',[])
        
    else
        ticklabels_new = cell(size(Labels));
        for i = 1:length(Labels)
            ticklabels_new{i} = ['\color{black} ' Labels{i}];
        end
        % set the tick labels
        set(gca, 'XTickLabel', ticklabels_new,'FontSize',20,'FontName', 'Arial');
    end

end
