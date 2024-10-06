%Function to plot bar plots and scatter for condition averages

function [] = avg_bar_scatter_simulations(results, sum_data, Labels, cols, col_idx, phase, ylim_lower, mod_num)
    y_ticks = ylim_lower:.2:1;

    %plot p_correct broken up for conditions
    if sum(results(:,8) == 1,1) > sum(results(:,8) == 2,1)
        data_to_plot = NaN(sum(results(:,8) == 1,1),4);
    else
        data_to_plot = NaN(sum(results(:,8) == 2,1),4);
    end
    
%sorted by condition
%     data_to_plot(1:sum(results(:,8) == 1),1) = sum_data(results(:,8) == 1,col_idx(1));
%     data_to_plot(1:sum(results(:,8) == 2),2) = sum_data(results(:,8) == 1,col_idx(2));
%     data_to_plot(1:sum(results(:,8) == 1),3) = sum_data(results(:,8) == 2,col_idx(1));
%     data_to_plot(1:sum(results(:,8) == 2),4) = sum_data(results(:,8) == 2,col_idx(2));

    %sorted by trial type
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
    b = bar(1:4,nanmean(data_to_plot),'FaceColor','flat', 'BarWidth', 0.6, 'linewidth', 2.5, 'FaceAlpha',.5); hold all;
%     if phase == 1
%         b.CData(1,:) = cols.y;
%         b.CData(2,:) = cols.ytrend;
%         b.CData(3,:) = cols.b;
%         b.CData(4,:) = cols.btrend;
%     else
%         b.CData(1,:) = cols.green;
%         b.CData(2,:) = cols.k;
%         b.CData(3,:) = cols.green;
%         b.CData(4,:) = cols.k;
%     end
    if phase == 1
        b.CData(1,:) = cols.y;
        b.CData(2,:) = cols.b;
        b.CData(3,:) = cols.ytrend;
        b.CData(4,:) = cols.btrend;
    else
        b.CData(1,:) = cols.green;
        b.CData(2,:) = cols.green;
        b.CData(3,:) = cols.k;
        b.CData(4,:) = cols.k;
    end
    
    set(findobj(gca,'type','line'),'linew',2)

    errorbar(positions(1),nanmean(data_to_plot(:,1)), ci95(1), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
    errorbar(positions(2),nanmean(data_to_plot(:,2)), ci95(2), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
    errorbar(positions(3),nanmean(data_to_plot(:,3)), ci95(3), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
    errorbar(positions(4),nanmean(data_to_plot(:,4)), ci95(4), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);     
    
    %ylim([0 1]);  
    ylim([ylim_lower 1]);
    %xlim([0 5]);
    xlim([0.5 4.5]);
    box off
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    ybounds = ylim;
    yline(0.5,'--','linewidth',2.5)
%     set(gca,'YTick',[0 .25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
    set(gca,'YTick',y_ticks, 'FontSize',26,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'xtick',1:4)

    set(gcf,'color','w');
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)
    
    if mod_num == 4
         
        set(gca,'XTickLabel', Labels, 'FontSize',26,'FontName', 'Arial');
        set(gca,'XTickLabelRotation', 45);
        ticklabels_new = cell(size(Labels));
        for i = 1:length(Labels)
            ticklabels_new{i} = ['\color{black} ' Labels{i}];
        end
        % set the tick labels
        set(gca, 'XTickLabel', ticklabels_new,'FontSize',26,'FontName', 'Arial');
 
    end
   
    % prepend a color for each tick label
    LabelsY = get(gca,'YTickLabel');
    ticklabels_ynew = cell(size(LabelsY));
    for i = 1:length(LabelsY)
        ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_ynew,'FontSize',20,'FontName', 'Arial');


end
