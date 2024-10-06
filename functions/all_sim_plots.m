%create difference variables
if mod_num == 1
    clear diffs

    diffs(:,1) = nanmean(exp_res_to_plot(results(:,8)==1,5)) - nanmean(exp_res_to_plot(results(:,8)==2,5));
    diffs(:,2) = nanmean(exp_res_to_plot(results(:,8)==1,6)) - nanmean(exp_res_to_plot(results(:,8)==2,6));
    diffs(:,3) = nanmean(inf_to_plot(results(:,8)==1,5)) - nanmean(inf_to_plot(results(:,8)==2,5));
    diffs(:,4) = nanmean(inf_to_plot(results(:,8)==1,6)) - nanmean(inf_to_plot(results(:,8)==2,6));

    for k = 1:2
        diffs(2,k) = nanstd(exp_res_to_plot(results(:,8)==1,4+k)) / sqrt(length(exp_res_to_plot(results(:,8)==1,4+k))-1) + ...
                     nanstd(exp_res_to_plot(results(:,8)==2,4+k)) / sqrt(length(exp_res_to_plot(results(:,8)==2,4+k))-1);
    end

    for k = 1:2
        diffs(2,k+2) = nanstd(inf_to_plot(results(:,8)==1,4+k)) / sqrt(length(inf_to_plot(results(:,8)==1,4+k))-1) + ...
                       nanstd(inf_to_plot(results(:,8)==2,4+k)) / sqrt(length(inf_to_plot(results(:,8)==2,4+k))-1);
    end

elseif mod_num == 2
    clear diffs

    diffs_tmp = readtable('pdat_pdiff.csv');
    diffs(:,1) = mean(diffs_tmp.pdiff_mean(2));
    diffs(:,2) = mean(diffs_tmp.pdiff_mean(1));
    diffs(:,3) = mean(diffs_tmp.pdiff_mean(4));
    diffs(:,4) = mean(diffs_tmp.pdiff_mean(3));
    
    diffs(2,1) = (diffs_tmp.p1_std(2) / sqrt(1000)) + (diffs_tmp.p2_std(2) / sqrt(1000));
    diffs(2,2) = (diffs_tmp.p1_std(1) / sqrt(1000)) + (diffs_tmp.p2_std(1) / sqrt(1000));
    diffs(2,3) = (diffs_tmp.p1_std(4) / sqrt(1000)) + (diffs_tmp.p2_std(4) / sqrt(1000));
    diffs(2,4) = (diffs_tmp.p1_std(3) / sqrt(1000)) + (diffs_tmp.p2_std(3) / sqrt(1000));

elseif mod_num == 3
    
    diffs = sim_data.diffs;

end


positions = [1 2 3 4 5];   
pos = [positions-.3; ...
       positions+.3];
size_vec = [1 2 3 4 5];


nexttile
b = bar(1:4,diffs(1,:),'FaceColor','flat', 'BarWidth', 0.6, 'linewidth', 2.5, 'FaceAlpha',.5); hold all;
b.CData(1,:) = cols.green;
b.CData(2,:) = cols.k;
b.CData(3,:) = cols.green;
b.CData(4,:) = cols.k;

set(findobj(gca,'type','line'),'linew',2)

errorbar(positions(1),nanmean(diffs(1,1)), diffs(2,1), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
errorbar(positions(2),nanmean(diffs(1,2)), diffs(2,2), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
errorbar(positions(3),nanmean(diffs(1,3)), diffs(2,3), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);
errorbar(positions(4),nanmean(diffs(1,4)), diffs(2,4), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 4);

%ylim([0 1]);  
ylim([-0.1 .15]);
%xlim([0 5]);
xlim([0.5 4.5]);
box off
set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
ybounds = ylim;
%     set(gca,'YTick',[0 .25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'YTick',-0.1:.1:.15, 'FontSize',26,'FontName', 'Arial');
set(gca,'TickDir','out')
set(gca,'xtick',1:4)

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
%ylabel('\Delta Probability Correct','FontSize',20,'FontName', 'Arial')
h=gca; 
h.XAxis.TickLength = [0 0]; 
set(gca, 'XTickLabel', [])




