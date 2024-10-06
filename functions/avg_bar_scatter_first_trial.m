%Function to plot bar plots and scatter for condition averages

function [] = avg_bar_scatter_first_trial(data, markersize, cols, probe_type, Labels)

    %prepare data for plotting - within condition difference across trial
    %types
    data_to_plot(:,1) = data(data(:,1) == 1,4) - data(data(:,1) == 1,2);
    data_to_plot(:,3) = data(data(:,1) == 1,5) - data(data(:,1) == 1,3);

    data_to_plot(:,2) = data(data(:,1) == 2,4) - data(data(:,1) == 2,2);
    data_to_plot(:,4) = data(data(:,1) == 2,5) - data(data(:,1) == 2,3);

    %condition difference of differences
    diff_to_plot(:,1) = nanmean(data_to_plot(:,1)) - nanmean(data_to_plot(:,2));
    diff_to_plot(:,2) = nanmean(data_to_plot(:,3)) - nanmean(data_to_plot(:,4));

    

    positions = [1 2 3 4 5 6 7];   
    pos = [positions-.3; ...
           positions+.3];
    size_vec = [1 2 3 4 5 6 7];

    % %prepare confidence interval
    ci = 0.95 ;
    alpha = 1-ci;

    n = size(data_to_plot,1);
    T_multiplier = tinv(1-alpha/2, n-1);

    ci95 = T_multiplier*nanstd(data_to_plot)/sqrt(n);

    %plot with two different y-axes     
    hold on

    b = bar(1:6,[nanmean(data_to_plot(:,1:2)) diff_to_plot(1) nanmean(data_to_plot(:,3:4)) diff_to_plot(2)],'FaceColor','flat', 'BarWidth', 0.6, 'linewidth', 2.5, 'FaceAlpha',.5); hold all;
    b.CData([1 2 4 5],:) = 'w';
    b.CData(3,:) = cols.green;
    b.CData(6,:) = cols.k;

    errorbar(positions(1),nanmean(data_to_plot(:,1)), ci95(1), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
    errorbar(positions(2),nanmean(data_to_plot(:,2)), ci95(2), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
    errorbar(positions(4),nanmean(data_to_plot(:,3)), ci95(3), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
    errorbar(positions(5),nanmean(data_to_plot(:,4)), ci95(4), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);     
    ylim([-.30 .30]);     
    set(gca,'YTick',[-.30:.15:.30], 'FontSize',20,'FontName', 'Arial');       
    % prepend a color for each tick label
    LabelsY = get(gca,'YTickLabel');
    ticklabels_ynew = cell(size(LabelsY));
    for i = 1:length(LabelsY)
        ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_ynew,'FontSize',20,'FontName', 'Arial');
    ylabel('\Delta Probability Correct','FontSize',20,'FontName', 'Arial')

    
    yyaxis right
    scatter(linspace(pos(1,1), pos(2,1),length(data_to_plot(:,1))), data_to_plot(:,1), markersize, 'o', 'MarkerFaceColor',  cols.green, 'MarkerEdgeColor', cols.green,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
    scatter(linspace(pos(1,2), pos(2,2), length(data_to_plot(:,2))), data_to_plot(:,2), markersize, 'o', 'MarkerFaceColor',  cols.green, 'MarkerEdgeColor', cols.green,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
    scatter(linspace(pos(1,4), pos(2,4), length(data_to_plot(:,3))), data_to_plot(:,3), markersize, 'o', 'MarkerFaceColor',  cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 
    scatter(linspace(pos(1,5), pos(2,5), length(data_to_plot(:,4))), data_to_plot(:,4), markersize, 'o', 'MarkerFaceColor',  cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',0.001,'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha',.5); 

    ylim([-1 1]);  

    hold off
    
    %xlim([0 5]);
    xlim([0.5 6.5]);
    box off
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)

    ybounds = ylim;
    set(gca,'TickDir','out')
    set(gca,'xtick',1:6)
    set(gca,'XTickLabel', Labels, 'FontSize',20,'FontName', 'Arial');
    set(gca,'XTickLabelRotation', 45);
    set(gcf,'color','w');
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)

    if probe_type == 2
%        set(gca,'YTick',[])
%        ax1 = gca;                  
%        ax1.YAxis.Visible = 'off';   % remove y-axis

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
