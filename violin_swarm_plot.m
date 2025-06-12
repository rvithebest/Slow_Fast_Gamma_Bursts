function violin_swarm_plot(median_gatherer_sg,median_gatherer_fg)
    data=[(median_gatherer_sg)',(median_gatherer_fg)'];
    group=categorical({'Slow \gamma','Fast \gamma'});
    group=repmat(group,size(data,1),1);
    violinplot(group,data)
    v=violinplot(group,data);
    blue_color=[0 0 1];
    orange_color=[1 0.5 0];
    v(1).FaceColor=blue_color;
    v(2).FaceColor=orange_color;
    ylabel('Burst lengths(s)')
    %set(gca,'FontSize',25)
    hold on;
    % swarmchart plot- slow gamma - blue, fast gamma - orange
    color=[blue_color;orange_color];
    swarmchart(group,data,[],color,'filled')
    hold on;
    med_median_gatherer_sg=median(median_gatherer_sg);
    sem_median_gatherer_sg=getSEMedian(median_gatherer_sg,1000);
    med_median_gatherer_fg=median(median_gatherer_fg);
    sem_median_gatherer_fg=getSEMedian(median_gatherer_fg,1000);
    % errorbar plot
    errorbar(2,med_median_gatherer_sg,sem_median_gatherer_sg,'g','LineWidth',2,'MarkerSize',15)
    errorbar(1,med_median_gatherer_fg,sem_median_gatherer_fg,'g','LineWidth',2,'MarkerSize',15)
    %%%%%%%%%%%%%%%%%% BOXPLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % superimpose a boxplot on the violin plot
    % data=[(median_gatherer_sg)',(median_gatherer_fg)'];
    % group=({'Slow gamma','Fast gamma'});
    % boxchart(data)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2-sample unpaired t-test
    [h,p,ci,stats]=ttest2(median_gatherer_sg,median_gatherer_fg,'VarType','unequal');
    % display in text box above the plot- significance- format **** (p-value)
    if p<0.0001
        %text(1.5,0.5,['***(',num2str(p),')'],'FontSize',15)
        significance='****';
    elseif p<0.001
        %text(1.5,0.5,['***(',num2str(p),')'],'FontSize',15)
        significance='***';
    elseif p<0.01
        %text(1.5,0.5,['**(',num2str(p),')'],'FontSize',15)
        significance='**';
    elseif p<0.05
        %text(1.5,0.5,['*(',num2str(p),')'],'FontSize',15)
        significance='*';
    else
        %text(1.5,0.5,['N.S.(',num2str(p),')'],'FontSize',15)
        significance='N.S.';
    end
   % Define x positions for the plots (e.g., positions 1 and 2)
    x_positions = [1, 2];
    % Calculate y positions for the line
    y_max = max([max(median_gatherer_sg), max(median_gatherer_fg)]);
    y_line = y_max + 0.05; % Add a buffer above the max value
    % Draw the significance line
    line(x_positions, [y_line y_line], 'Color', 'k', 'LineWidth', 2);
    text(mean(x_positions), y_line + 0.005, significance, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 14, 'FontWeight', 'bold');
end