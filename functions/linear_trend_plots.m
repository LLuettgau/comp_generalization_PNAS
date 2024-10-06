
%run some checks
nanmean(all_corr_resp_prim(results(:,8) == 1,2:3))
nanmean(all_corr_resp_prim(results(:,8) == 2,2:3))
nanmean(all_corr_resp_sec(results(:,8) == 1,2:3))
nanmean(all_corr_resp_sec(results(:,8) == 2,2:3))

nanmean(inf_all_corr_resp_prim(results(:,8) == 1,2:3))
nanmean(inf_all_corr_resp_prim(results(:,8) == 2,2:3))
nanmean(inf_all_corr_resp_sec(results(:,8) == 1,2:3))
nanmean(inf_all_corr_resp_sec(results(:,8) == 2,2:3))


nanmean(trans_corr_resp_prim(results(:,8) == 1,2:3))
nanmean(trans_corr_resp_prim(results(:,8) == 2,2:3))
nanmean(trans_corr_resp_sec(results(:,8) == 1,2:3))
nanmean(trans_corr_resp_sec(results(:,8) == 2,2:3))


nanmean(trans_inf_corr_resp_prim(results(:,8) == 1,2:3))
nanmean(trans_inf_corr_resp_prim(results(:,8) == 2,2:3))
nanmean(trans_inf_corr_resp_sec(results(:,8) == 1,2:3))
nanmean(trans_inf_corr_resp_sec(results(:,8) == 2,2:3))

%plot binned data
h = figure;
t = tiledlayout(2,1,'TileSpacing','loose', 'Padding', 'tight');

num_bins = 4;

nexttile
all_trend_data_prior = plot_bins(corr_resp_all(:,2:end),inf_corr_resp_all(:,2:end),num_bins,[cols.dred; cols.lyellow]);
title('Prior Learning','FontSize',35,'FontName', 'Arial')


nexttile
all_trend_data_trans = plot_bins(trans_corr_resp_all(:,2:end),trans_inf_corr_resp_all(:,2:end),num_bins,[cols.medorange; cols.cream]);
title('Transfer Learning','FontSize',35,'FontName', 'Arial')
xlabel('Binned Trials', 'FontSize',30, ...
    'Color','k')
ylabel(t, 'Probability Correct','FontSize',30,'FontName', 'Arial')

% %save linear trend data to csv file (for analysis in R)
% all_trend_data = [repmat(1:4,1,4)' [ones(4,1); ones(4,1)+1'; ones(4,1)+2'; ones(4,1)+3'] [all_trend_data_prior'; all_trend_data_trans']];
% cd(save_dir)
% 
% varNames_lin_trend = {'bin','tp_pt', 'correct'};
% 
% lin_trend_data_table = array2table(all_trend_data,'VariableNames',varNames_lin_trend);
% writetable(lin_trend_data_table,['lin_trend_data',num2str(sample),'.csv'],'Delimiter',',');


cd(plotpath)


if save_plots == 1
    
    screen_size = get(0, 'ScreenSize');
    origSize = get(h, 'Position'); % grab original on screen size
    set(h, 'Position', [0 0 screen_size(3)/2 screen_size(4) ] ); %set to scren size
    set(h,'PaperPositionMode','auto') %set paper pos for printing
    % saveas(h, 'fig_trajectories.png') % save figure
    ax = h;
    exportgraphics(ax,'fig_trajectories_binned_paper.png','Resolution',1200)
    set(h,'Position', origSize) %set back to original dimensions
    
end








%linear trend plots

% Fit linear trends to data
%overall exp - prior
linearCoefficients_exp = polyfit(1:trial_num_prior,nanmean(corr_resp_all(:,2:end)),1); % First degree linear fit
yFit_exp = polyval(linearCoefficients_exp, 1:trial_num_prior);

%overall inf - prior
linearCoefficients_inf = polyfit(1:trial_num_prior,nanmean(inf_corr_resp_all(:,2:end)),1); % First degree linear fit
yFit_inf = polyval(linearCoefficients_inf, 1:trial_num_prior);


%overall exp - transfer
linearCoefficients_exp_trans = polyfit(1:trial_num_trans,nanmean(trans_corr_resp_all(:,2:end)),1); % First degree linear fit
yFit_exp_trans = polyval(linearCoefficients_exp_trans, 1:trial_num_trans);

%overall inf - transfer
linearCoefficients_inf_trans = polyfit(1:trial_num_trans,nanmean(trans_inf_corr_resp_all(:,2:end)),1); % First degree linear fit
yFit_inf_trans = polyval(linearCoefficients_inf_trans, 1:trial_num_trans);


smooth_param = 6;

h = figure;
t = tiledlayout(2,1,'TileSpacing','loose', 'Padding', 'tight');

nexttile
stdshade(corr_resp_all(:,2:size(corr_resp_all,2)), .4, cols.dred, [], smooth_param); hold on;
stdshade(inf_corr_resp_all(:,2:size(inf_corr_resp_all,2)), .4, cols.lyellow, [], smooth_param);
set(findobj(gca,'type','line'),'linew',2);
p1 = plot(1:trial_num_prior,yFit_exp, '--','Color', cols.dgrey);
p2 = plot(1:trial_num_prior,yFit_inf, '--','Color', cols.dgrey);
p1.Color(4) = 0.5;
p2.Color(4) = 0.5;
hold off
box off
xbounds = [0 size(corr_resp_all,2)-1];
ylim([0.5 .85]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.5:.10:.80], 'FontSize',25,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):smooth_param:xbounds(2), 'FontSize',25,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
title('Prior Learning','FontSize',35,'FontName', 'Arial')

xtickangle(45)

nexttile
stdshade(trans_corr_resp_all(:,2:size(trans_corr_resp_all,2)), .4, cols.medorange, [], smooth_param); hold on;
stdshade(trans_inf_corr_resp_all(:,2:size(trans_inf_corr_resp_all,2)), .4, cols.cream, [], smooth_param);
set(findobj(gca,'type','line'),'linew',2);
p1 = plot(1:trial_num_trans,yFit_exp_trans, '--','Color', cols.dgrey);
p2 = plot(1:trial_num_trans,yFit_inf_trans, '--','Color', cols.dgrey);
p1.Color(4) = 0.5;
p2.Color(4) = 0.5;
hold off
box off
xbounds = [0 size(trans_corr_resp_all,2)-1];
ylim([0.5 .85]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.5:.10:.80], 'FontSize',25,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):smooth_param:xbounds(2), 'FontSize',25,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
xlabel('Trial Number', 'FontSize',30, ...
    'Color','k')
xtickangle(45)
title('Transfer Learning','FontSize',35,'FontName', 'Arial')
ylabel(t, 'Probability Correct','FontSize',30,'FontName', 'Arial')





% if save_plots == 1
%     
%     screen_size = get(0, 'ScreenSize');
%     origSize = get(h, 'Position'); % grab original on screen size
%     set(h, 'Position', [0 0 screen_size(3)/2 screen_size(4) ] ); %set to scren size
%     set(h,'PaperPositionMode','auto') %set paper pos for printing
%     % saveas(h, 'fig_trajectories.png') % save figure
%     ax = h;
%     exportgraphics(ax,'fig_trajectories_paper.png','Resolution',1200)
%     set(h,'Position', origSize) %set back to original dimensions
%     
% end










% Fit linear trends to data
%overall exp
linearCoefficients_exp = polyfit(1:trial_num_prior,nanmean(corr_resp_all(:,2:end)),1); % First degree linear fit
yFit_exp = polyval(linearCoefficients_exp, 1:trial_num_prior);

%overall inf
linearCoefficients_inf = polyfit(1:trial_num_prior,nanmean(inf_corr_resp_all(:,2:end)),1); % First degree linear fit
yFit_inf = polyval(linearCoefficients_inf, 1:trial_num_prior);

h = figure;
%subplot(1,3,2)
subplot(1,2,1)
stdshade(corr_resp_all(:,2:size(corr_resp_all,2)), .4, cols.lblue, [], 5); hold on;
stdshade(inf_corr_resp_all(:,2:size(inf_corr_resp_all,2)), .4, cols.green, [], 5);
set(findobj(gca,'type','line'),'linew',2);
plot(1:trial_num_prior,yFit_exp, '--','Color', cols.btrend)
plot(1:trial_num_prior,yFit_inf, '--','Color', cols.gtrend)
hold off
box off
xbounds = [0 size(corr_resp_all,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):10:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
ylabel('Probability of Correct Answers', 'FontSize',20,...
    'Color','k')
title('Average correct responses - Prior')
xtickangle(45)
legend('','Seq. Probes', '','Inf Probes','Location','NorthWestOutside', 'FontSize',12)


% Fit linear trends to data - exp
%overall primary
linearCoefficients_prim_exp = polyfit(1:trial_num_prior/2,nanmean(corr_resp_prim(:,2:end)),1); % First degree linear fit
yFit_prim_exp = polyval(linearCoefficients_prim_exp, 1:trial_num_prior/2);

%overall secondary
linearCoefficients_sec_exp = polyfit(1:trial_num_prior/2,nanmean(corr_resp_sec(:,2:end)),1); % First degree linear fit
yFit_sec_exp = polyval(linearCoefficients_sec_exp, 1:trial_num_prior/2);

%%average p_corr for trial types
subplot(2,2,2)
stdshade(corr_resp_prim(:,2:size(corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
stdshade(corr_resp_sec(:,2:size(corr_resp_sec,2)), .1, cols.b, [], 1);
set(findobj(gca,'type','line'),'linew',2);
plot(1:trial_num_prior/2,yFit_prim_exp, '--','Color', cols.ytrend)
plot(1:trial_num_prior/2,yFit_sec_exp, '--','Color', cols.btrend)
hold off
box off
xbounds = [0 size(corr_resp_prim,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
% ylabel('Probability of Correct Answers', 'FontSize',20,...
%     'Color','k')
title('Seq Learning: Trial type')
legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
xtickangle(45)


% Fit linear trends to data - inference
%overall primary
linearCoefficients_prim = polyfit(1:trial_num_prior/2,nanmean(inf_corr_resp_prim(:,2:end)),1); % First degree linear fit
yFit_prim = polyval(linearCoefficients_prim, 1:trial_num_prior/2);

%overall secondary
linearCoefficients_sec = polyfit(1:trial_num_prior/2,nanmean(inf_corr_resp_sec(:,2:end)),1); % First degree linear fit
yFit_sec = polyval(linearCoefficients_sec, 1:trial_num_prior/2);

%%average p_corr for trial types - inference
subplot(2,2,4)
stdshade(inf_corr_resp_prim(:,2:size(inf_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
stdshade(inf_corr_resp_sec(:,2:size(inf_corr_resp_sec,2)), .1, cols.b, [], 1);
set(findobj(gca,'type','line'),'linew',2);
plot(1:trial_num_prior/2,yFit_prim, '--','Color', cols.ytrend)
plot(1:trial_num_prior/2,yFit_sec, '--','Color', cols.btrend)
hold off
box off
xbounds = [0 size(inf_corr_resp_prim,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
% ylabel('Probability of Correct Answers', 'FontSize',20,...
%     'Color','k')
title('Inference: Trial type')
legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
xtickangle(45)

if save_plots == 1
    
    screen_size = get(0, 'ScreenSize');
    origSize = get(h, 'Position'); % grab original on screen size
    set(h, 'Position', [0 0 screen_size(3) screen_size(4)/2 ] ); %set to scren size
    set(h,'PaperPositionMode','auto') %set paper pos for printing
    % saveas(h, 'fig_trajectories.png') % save figure
    ax = h;
    exportgraphics(ax,'fig_trajectories.png','Resolution',1200)
    set(h,'Position', origSize) %set back to original dimensions
    
end

%% Trajectories in prior learning - primary/secondary - split up for conditions

% Fit linear trends to data - exp
%overall primary
linearCoefficients_prim_exp_prior_cond1 = polyfit(1:trial_num_prior/2,nanmean(corr_resp_prim(results_unsorted(:,8) == 1,2:end)),1); % First degree linear fit
yFit_prim_exp_prior_cond1 = polyval(linearCoefficients_prim_exp_prior_cond1, 1:trial_num_prior/2);

linearCoefficients_prim_exp_prior_cond2 = polyfit(1:trial_num_prior/2,nanmean(corr_resp_prim(results_unsorted(:,8) == 2,2:end)),1); % First degree linear fit
yFit_prim_exp_prior_cond2 = polyval(linearCoefficients_prim_exp_prior_cond2, 1:trial_num_prior/2);

%overall secondary
linearCoefficients_sec_exp_prior_cond1 = polyfit(1:trial_num_prior/2,nanmean(corr_resp_sec(results_unsorted(:,8) == 1,2:end)),1); % First degree linear fit
yFit_sec_exp_prior_cond1 = polyval(linearCoefficients_sec_exp_prior_cond1, 1:trial_num_prior/2);

linearCoefficients_sec_exp_prior_cond2 = polyfit(1:trial_num_prior/2,nanmean(corr_resp_sec(results_unsorted(:,8) == 2,2:end)),1); % First degree linear fit
yFit_sec_exp_prior_cond2 = polyval(linearCoefficients_sec_exp_prior_cond2, 1:trial_num_prior/2);

%%average p_corr for trial types
%condition 1
h = figure;
subplot(2,2,1)
stdshade(corr_resp_prim(results_unsorted(:,8) == 1,2:size(corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
stdshade(corr_resp_sec(results_unsorted(:,8) == 1,2:size(corr_resp_sec,2)), .1, cols.b, [], 1);
set(findobj(gca,'type','line'),'linew',2);
plot(1:trial_num_prior/2,yFit_prim_exp_prior_cond1, '--','Color', cols.ytrend)
plot(1:trial_num_prior/2,yFit_sec_exp_prior_cond1, '--','Color', cols.btrend)
hold off
box off
xbounds = [0 size(corr_resp_prim,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
% ylabel('Probability of Correct Answers', 'FontSize',20,...
%     'Color','k')
title('Seq Learning: Trial type - Prior, cond1')
legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
xtickangle(45)


%%average p_corr for trial types
%condition 2
subplot(2,2,3)
stdshade(corr_resp_prim(results_unsorted(:,8) == 2,2:size(corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
stdshade(corr_resp_sec(results_unsorted(:,8) == 2,2:size(corr_resp_sec,2)), .1, cols.b, [], 1);
set(findobj(gca,'type','line'),'linew',2);
plot(1:trial_num_prior/2,yFit_prim_exp_prior_cond2, '--','Color', cols.ytrend)
plot(1:trial_num_prior/2,yFit_sec_exp_prior_cond2, '--','Color', cols.btrend)
hold off
box off
xbounds = [0 size(corr_resp_prim,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
% ylabel('Probability of Correct Answers', 'FontSize',20,...
%     'Color','k')
title('Seq Learning: Trial type - Prior, cond2')
legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
xtickangle(45)



%%average p_corr for trial types - inference, transfer
%overall primary
linearCoefficients_prim_inf_prior_cond1 = polyfit(1:trial_num_prior/2,nanmean(inf_corr_resp_prim(results_unsorted(:,8) == 1,2:end)),1); % First degree linear fit
yFit_prim_inf_prior_cond1 = polyval(linearCoefficients_prim_inf_prior_cond1, 1:trial_num_prior/2);

linearCoefficients_prim_inf_prior_cond2 = polyfit(1:trial_num_prior/2,nanmean(inf_corr_resp_prim(results_unsorted(:,8) == 2,2:end)),1); % First degree linear fit
yFit_prim_inf_prior_cond2 = polyval(linearCoefficients_prim_inf_prior_cond2, 1:trial_num_prior/2);

%overall secondary
linearCoefficients_sec_trans_inf_cond1 = polyfit(1:trial_num_prior/2,nanmean(inf_corr_resp_sec(results_unsorted(:,8) == 1,2:end)),1); % First degree linear fit
yFit_sec_inf_prior_cond1 = polyval(linearCoefficients_sec_trans_inf_cond1, 1:trial_num_prior/2);

linearCoefficients_sec_inf_prior_cond2 = polyfit(1:trial_num_prior/2,nanmean(inf_corr_resp_sec(results_unsorted(:,8) == 2,2:end)),1); % First degree linear fit
yFit_sec_inf_prior_cond2 = polyval(linearCoefficients_sec_inf_prior_cond2, 1:trial_num_prior/2);

%condition 1
subplot(2,2,2)
stdshade(inf_corr_resp_prim(results_unsorted(:,8) == 1,2:size(inf_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
stdshade(inf_corr_resp_sec(results_unsorted(:,8) == 1,2:size(inf_corr_resp_sec,2)), .1, cols.b, [], 1);
set(findobj(gca,'type','line'),'linew',2);
plot(1:trial_num_prior/2,yFit_prim_inf_prior_cond1, '--','Color', cols.ytrend)
plot(1:trial_num_prior/2,yFit_sec_inf_prior_cond1, '--','Color', cols.btrend)
hold off
box off
xbounds = [0 size(inf_corr_resp_prim,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
% ylabel('Probability of Correct Answers', 'FontSize',20,...
%     'Color','k')
title('Inference: Trial type - Prior, cond1')
legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
xtickangle(45)

%condition 2
subplot(2,2,4)
stdshade(inf_corr_resp_prim(results_unsorted(:,8) == 2,2:size(inf_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
stdshade(inf_corr_resp_sec(results_unsorted(:,8) == 2,2:size(inf_corr_resp_sec,2)), .1, cols.b, [], 1);
set(findobj(gca,'type','line'),'linew',2);
plot(1:trial_num_prior/2,yFit_prim_inf_prior_cond2, '--','Color', cols.ytrend)
plot(1:trial_num_prior/2,yFit_sec_inf_prior_cond2, '--','Color', cols.btrend)
hold off
box off
xbounds = [0 size(inf_corr_resp_prim,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
% ylabel('Probability of Correct Answers', 'FontSize',20,...
%     'Color','k')
title('Inference: Trial type - Prior, cond2')
legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
xtickangle(45)


if save_plots == 1
    
    screen_size = get(0, 'ScreenSize');
    origSize = get(h, 'Position'); % grab original on screen size
    set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size
    set(h,'PaperPositionMode','auto') %set paper pos for printing
    % saveas(h, 'fig_trajectories.png') % save figure
    ax = h;
    exportgraphics(ax,'fig_trajectories_trial_type_prior.png','Resolution',1200)
    set(h,'Position', origSize) %set back to original dimensions
    
end




%% Trajectories in transfer overall - split up for conditions

%average p_corr - exp and inference probes - transfer learning

% Fit linear trends to data
%overall exp
linearCoefficients_trans_exp_cond1 = polyfit(1:36,nanmean(trans_corr_resp_all(results_unsorted(:,8) == 1,2:end)),1); % First degree linear fit
yFit_trans_exp_cond1 = polyval(linearCoefficients_trans_exp_cond1, 1:36);

linearCoefficients_trans_exp_cond2 = polyfit(1:36,nanmean(trans_corr_resp_all(results_unsorted(:,8) == 2,2:end)),1); % First degree linear fit
yFit_trans_exp_cond2 = polyval(linearCoefficients_trans_exp_cond2, 1:36);

%overall inf
linearCoefficients_trans_inf_cond1 = polyfit(1:36,nanmean(trans_inf_corr_resp_all(results_unsorted(:,8) == 1,2:end)),1); % First degree linear fit
yFit_trans_inf_cond1 = polyval(linearCoefficients_trans_inf_cond1, 1:36);

linearCoefficients_trans_inf_cond2 = polyfit(1:36,nanmean(trans_inf_corr_resp_all(results_unsorted(:,8) == 2,2:end)),1); % First degree linear fit
yFit_trans_inf_cond2 = polyval(linearCoefficients_trans_inf_cond2, 1:36);


h = figure;
%subplot(1,3,2)
subplot(1,2,1)
stdshade(trans_corr_resp_all(results_unsorted(:,8) == 1,2:size(trans_corr_resp_all,2)), .4, cols.lblue, [], 5); hold on;
stdshade(trans_inf_corr_resp_all(results_unsorted(:,8) == 1,2:size(trans_inf_corr_resp_all,2)), .4, cols.green, [], 5);
set(findobj(gca,'type','line'),'linew',2);
plot(1:36,yFit_trans_exp_cond1, '--','Color', cols.btrend)
plot(1:36,yFit_trans_inf_cond1, '--','Color', cols.gtrend)
hold off
box off
xbounds = [0 size(trans_corr_resp_all,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
ylabel('Probability of Correct Answers', 'FontSize',20,...
    'Color','k')
title('Average correct responses - Transfer, cond1')
xtickangle(45)
%legend('','Seq. Probes', '','Inf Probes','Location','NorthWestOutside', 'FontSize',12)

subplot(1,2,2)
stdshade(trans_corr_resp_all(results_unsorted(:,8) == 2,2:size(trans_corr_resp_all,2)), .4, cols.lblue, [], 5); hold on;
stdshade(trans_inf_corr_resp_all(results_unsorted(:,8) == 2,2:size(trans_inf_corr_resp_all,2)), .4, cols.green, [], 5);
set(findobj(gca,'type','line'),'linew',2);
plot(1:36,yFit_trans_exp_cond2, '--','Color', cols.btrend)
plot(1:36,yFit_trans_inf_cond2, '--','Color', cols.gtrend)
hold off
box off
xbounds = [0 size(trans_corr_resp_all,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
ylabel('Probability of Correct Answers', 'FontSize',20,...
    'Color','k')
title('Average correct responses - Transfer, cond2')
xtickangle(45)
%legend('','Seq. Probes', '','Inf Probes','Location','NorthWestOutside', 'FontSize',12)

if save_plots == 1
    
    screen_size = get(0, 'ScreenSize');
    origSize = get(h, 'Position'); % grab original on screen size
    set(h, 'Position', [0 0 screen_size(3) screen_size(4)/2 ] ); %set to scren size
    set(h,'PaperPositionMode','auto') %set paper pos for printing
    % saveas(h, 'fig_trajectories.png') % save figure
    ax = h;
    exportgraphics(ax,'fig_trajectories_transfer.png','Resolution',1200)
    set(h,'Position', origSize) %set back to original dimensions
end


%% Trajectories in transfer - primary/secondary - split up for conditions
% Fit linear trends to data - exp
%overall primary
linearCoefficients_prim_exp_trans_cond1 = polyfit(1:18,nanmean(trans_corr_resp_prim(results_unsorted(:,8) == 1,2:end)),1); % First degree linear fit
yFit_prim_exp_trans_cond1 = polyval(linearCoefficients_prim_exp_trans_cond1, 1:18);

[linearCoefficients_prim_exp_trans_cond2, S] = polyfit(1:18,nanmean(trans_corr_resp_prim(results_unsorted(:,8) == 2,2:end)),1); % First degree linear fit
yFit_prim_exp_trans_cond2 = polyval(linearCoefficients_prim_exp_trans_cond2, 1:18);
CI = polyparci(linearCoefficients_prim_exp_trans_cond2,S);

%overall secondary
linearCoefficients_sec_exp_trans_cond1 = polyfit(1:18,nanmean(trans_corr_resp_sec(results_unsorted(:,8) == 1,2:end)),1); % First degree linear fit
yFit_sec_exp_trans_cond1 = polyval(linearCoefficients_sec_exp_trans_cond1, 1:18);

linearCoefficients_sec_exp_trans_cond2 = polyfit(1:18,nanmean(trans_corr_resp_sec(results_unsorted(:,8) == 2,2:end)),1); % First degree linear fit
yFit_sec_exp_trans_cond2 = polyval(linearCoefficients_sec_exp_trans_cond2, 1:18);

%%average p_corr for trial types
%condition 1
h = figure;
subplot(2,2,1)
stdshade(trans_corr_resp_prim(results_unsorted(:,8) == 1,2:size(trans_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
stdshade(trans_corr_resp_sec(results_unsorted(:,8) == 1,2:size(trans_corr_resp_sec,2)), .1, cols.b, [], 1);
set(findobj(gca,'type','line'),'linew',2);
plot(1:18,yFit_prim_exp_trans_cond1, '--','Color', cols.ytrend)
plot(1:18,yFit_sec_exp_trans_cond1, '--','Color', cols.btrend)
hold off
box off
xbounds = [0 size(trans_corr_resp_prim,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
% ylabel('Probability of Correct Answers', 'FontSize',20,...
%     'Color','k')
title('Seq Learning: Trial type - Transfer, cond1')
legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
xtickangle(45)


%%average p_corr for trial types
%condition 2
subplot(2,2,3)
stdshade(trans_corr_resp_prim(results_unsorted(:,8) == 2,2:size(trans_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
stdshade(trans_corr_resp_sec(results_unsorted(:,8) == 2,2:size(trans_corr_resp_sec,2)), .1, cols.b, [], 1);
set(findobj(gca,'type','line'),'linew',2);
plot(1:18,yFit_prim_exp_trans_cond2, '--','Color', cols.ytrend)
plot(1:18,yFit_sec_exp_trans_cond2, '--','Color', cols.btrend)
hold off
box off
xbounds = [0 size(trans_corr_resp_prim,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
% ylabel('Probability of Correct Answers', 'FontSize',20,...
%     'Color','k')
title('Seq Learning: Trial type - Transfer, cond2')
legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
xtickangle(45)



%%average p_corr for trial types - inference, transfer
%overall primary
linearCoefficients_prim_inf_trans_cond1 = polyfit(1:18,nanmean(trans_inf_corr_resp_prim(results_unsorted(:,8) == 1,2:end)),1); % First degree linear fit
yFit_prim_inf_trans_cond1 = polyval(linearCoefficients_prim_inf_trans_cond1, 1:18);

linearCoefficients_prim_inf_trans_cond2 = polyfit(1:18,nanmean(trans_inf_corr_resp_prim(results_unsorted(:,8) == 2,2:end)),1); % First degree linear fit
yFit_prim_inf_trans_cond2 = polyval(linearCoefficients_prim_inf_trans_cond2, 1:18);

%overall secondary
linearCoefficients_sec_inf_trans_cond1 = polyfit(1:18,nanmean(trans_inf_corr_resp_sec(results_unsorted(:,8) == 1,2:end)),1); % First degree linear fit
yFit_sec_inf_trans_cond1 = polyval(linearCoefficients_sec_inf_trans_cond1, 1:18);

linearCoefficients_sec_inf_trans_cond2 = polyfit(1:18,nanmean(trans_inf_corr_resp_sec(results_unsorted(:,8) == 2,2:end)),1); % First degree linear fit
yFit_sec_inf_trans_cond2 = polyval(linearCoefficients_sec_inf_trans_cond2, 1:18);

%condition 1
subplot(2,2,2)
stdshade(trans_inf_corr_resp_prim(results_unsorted(:,8) == 1,2:size(trans_inf_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
stdshade(trans_inf_corr_resp_sec(results_unsorted(:,8) == 1,2:size(trans_inf_corr_resp_sec,2)), .1, cols.b, [], 1);
set(findobj(gca,'type','line'),'linew',2);
plot(1:18,yFit_prim_inf_trans_cond1, '--','Color', cols.ytrend)
plot(1:18,yFit_sec_inf_trans_cond1, '--','Color', cols.btrend)
hold off
box off
xbounds = [0 size(trans_inf_corr_resp_prim,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
% ylabel('Probability of Correct Answers', 'FontSize',20,...
%     'Color','k')
title('Inference: Trial type - Transfer, cond1')
legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
xtickangle(45)

%condition 2
subplot(2,2,4)
stdshade(trans_inf_corr_resp_prim(results_unsorted(:,8) == 2,2:size(trans_inf_corr_resp_prim,2)), .1, cols.y, [], 1); hold on;
stdshade(trans_inf_corr_resp_sec(results_unsorted(:,8) == 2,2:size(trans_inf_corr_resp_sec,2)), .1, cols.b, [], 1);
set(findobj(gca,'type','line'),'linew',2);
plot(1:18,yFit_prim_inf_trans_cond2, '--','Color', cols.ytrend)
plot(1:18,yFit_sec_inf_trans_cond2, '--','Color', cols.btrend)
hold off
box off
xbounds = [0 size(trans_inf_corr_resp_prim,2)-1];
ylim([0.25 1]);
xlim(xbounds);
set(findobj(gca,'type','line'),'linew',3)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
set(gca,'YTick',[.25 .5 .75 1], 'FontSize',20,'FontName', 'Arial');
set(gca,'XTick',xbounds(1):5:xbounds(2), 'FontSize',20,'FontName', 'Arial');
yline(0.5,'k--');
set(gcf,'color','w');
set(gca,'TickDir','out')
% ylabel('Probability of Correct Answers', 'FontSize',20,...
%     'Color','k')
title('Inference: Trial type - Transfer, cond2')
legend('','Primary', '','Secondary','Location','NorthEastOutside', 'FontSize',10)
xtickangle(45)


if save_plots == 1
    
    screen_size = get(0, 'ScreenSize');
    origSize = get(h, 'Position'); % grab original on screen size
    set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size
    set(h,'PaperPositionMode','auto') %set paper pos for printing
    % saveas(h, 'fig_trajectories.png') % save figure
    ax = h;
    exportgraphics(ax,'fig_trajectories_trial_type_transfer.png','Resolution',1200)
    set(h,'Position', origSize) %set back to original dimensions
end

