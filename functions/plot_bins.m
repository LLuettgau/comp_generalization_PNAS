function all_trend_data = plot_bins(timeseries_data_exp,timeseries_data_inf,bins,cols)
        
    num_bins = bins;
    bin_start_idx = 1:size(timeseries_data_exp,2)/num_bins:size(timeseries_data_exp,2);
    bin_end_idx = size(timeseries_data_exp,2)/num_bins:size(timeseries_data_exp,2)/num_bins:size(timeseries_data_exp,2);
    Labels = [];
    for k = 1:num_bins
        Labels_tmp = {[num2str(bin_start_idx(k)) ' - ' num2str(bin_end_idx(k))]};
        Labels = [Labels_tmp; Labels];
    end
    Labels = flip(Labels);

    %exp probes
    binned_data_exp = zeros(num_bins, 1);
    for i = 1:num_bins
        binned_data_exp(i) = mean(nanmean(timeseries_data_exp(:,bin_start_idx(i):bin_end_idx(i))));
        error_data_exp(i) = mean(nanstd(timeseries_data_exp(:,bin_start_idx(i):bin_end_idx(i))));
    end
    
    %inf probes
    binned_data_inf = zeros(num_bins, 1);
    for i = 1:num_bins
        binned_data_inf(i) = mean(nanmean(timeseries_data_inf(:,bin_start_idx(i):bin_end_idx(i))));
        error_data_inf(i) = mean(nanstd(timeseries_data_inf(:,bin_start_idx(i):bin_end_idx(i))));
    end

    
    %prepare SEM
    n = size(timeseries_data_exp,1);
    error_bars_exp = error_data_exp/sqrt(n);
    error_bars_inf = error_data_inf/sqrt(n);
    

    positions = 1:num_bins+1;   

    
    
    % Mean values with error bars
    for i = 1:num_bins
        errorbar(positions(i),binned_data_exp(i), error_bars_exp(i), 'o', 'MarkerSize', 20, ...
            'MarkerEdgeColor',cols(1,:),'LineWidth',0.001,'MarkerFaceColor',cols(1,:), ...
            'color', cols(1,:), 'CapSize', 25, 'LineWidth', 3); hold on
        
        errorbar(positions(i),binned_data_inf(i), error_bars_inf(i), 'o', 'MarkerSize', 20, ...
            'MarkerEdgeColor',cols(2,:),'LineWidth',0.001,'MarkerFaceColor',cols(2,:), ...
            'color', cols(2,:), 'CapSize', 25, 'LineWidth', 3); hold on
    end
    %linear trend
    [linearCoefficients_exp, S_exp] = polyfit(1:length(binned_data_exp),binned_data_exp,1); % First degree linear fit
    yFit_exp = polyval(linearCoefficients_exp, 1:length(binned_data_exp));
   
    
    disp('experience probes linear trend')
    CI_exp = polyparci(linearCoefficients_exp,S_exp);
    linearCoefficients_exp(:,1)
    CI_exp(:,1)
    
    
    
    [linearCoefficients_inf, S_inf] = polyfit(1:length(binned_data_inf),binned_data_inf,1); % First degree linear fit
    yFit_inf = polyval(linearCoefficients_inf, 1:length(binned_data_inf));
    
    disp('inference probes linear trend')
    CI_inf = polyparci(linearCoefficients_inf,S_inf);
    linearCoefficients_inf(:,1)
    CI_inf(:,1)
    
    

    p1 = plot(1:length(binned_data_exp),yFit_exp, '--','Color', [0.2,0.2,0.21]);
    p2 = plot(1:length(binned_data_inf),yFit_inf, '--','Color', [0.784,0.784,0.784]);
    p1.Color(4) = 0.5;
    p2.Color(4) = 0.5;
    set(findobj(gca,'type','line'),'linew',3)

    hold off
    ylim([0.5 .85]);
    xlim([0.5 num_bins + .5]);
    box off
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    set(gca,'YTick',[.5:.10:.80], 'FontSize',25,'FontName', 'Arial');

    ybounds = ylim;
    %yline(0.5,'--')
    %set(gca,'YTick',[0 .25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'xtick',1:num_bins)
    set(gca,'XTickLabel', Labels, 'FontSize',20,'FontName', 'Arial');
    set(gca,'XTickLabelRotation', 45);
    set(gcf,'color','w');
    set(gca,'ycolor',[0 0 0])
    set(gca,'xcolor',[0 0 0])

    all_trend_data = [binned_data_exp' binned_data_inf'];
end
