function violin_swarm_plot_paired(median_gatherer_sg,median_gatherer_fg,axis_ori)
    data=[(median_gatherer_sg)',(median_gatherer_fg)'];
    len=length(median_gatherer_sg);
    group=categorical({'Slow \gamma','Fast \gamma'});
    group=repmat(group,size(data,1),1);
    violinplot(group,data)
    v=violinplot(group,data);
    blue_color=[0 0 1];
    orange_color=[1 0.5 0];
    v(1).FaceColor=blue_color;
    v(2).FaceColor=orange_color;
    if axis_ori
       ylabel('Median Time center(s)');
    else
       ylabel('Median Burst length(s)');
    end
    %set(gca,'FontSize',25)
    hold on;
    % swarmchart plot- slow gamma - blue, fast gamma - orange
    color=[repmat(blue_color,len,1);repmat(orange_color,len,1)];
    % get coordinates of x axis
    x=[2*ones(1,len),ones(1,len)];
    data=[(median_gatherer_sg),(median_gatherer_fg)];
    s=swarmchart(x,data,'filled');
    s_struct=struct(s);
    swarmchart(x,data,[],color,'filled');
    a1=s_struct.XYZJittered(1:len,[1,2]);
    a2=s_struct.XYZJittered(1+len:end,[1,2]);
    for ii=1:length(a1)
       hold on;
       if a2(ii,2) > a1(ii,2)
            lineColor =  [1.00 0.00 1.00 0.3];  % Magenta
       else
            lineColor =  [0.35, 0, 0.6, 0.3];  % Dark violet
       end
       line([a1(ii,1), a2(ii,1)], [a1(ii,2), a2(ii,2)], 'Color', lineColor, 'LineWidth', 1)
    end
    hold on;
    if axis_ori==2
        med_median_gatherer_sg=mean(median_gatherer_sg);
        sem_median_gatherer_sg=std(median_gatherer_sg)/sqrt(length(median_gatherer_sg));
        med_median_gatherer_fg=mean(median_gatherer_fg);
        sem_median_gatherer_fg=std(median_gatherer_fg)/sqrt(length(median_gatherer_fg));
        % 2-sample paired t-tests
        [h,p,ci,stats]=ttest(median_gatherer_sg,median_gatherer_fg,'Alpha',0.05,'Tail','right');
        disp(["p-value equals to",string(p)])
        ylabel("Mean Onset time (s)")
    else
        med_median_gatherer_sg=median(median_gatherer_sg);
        sem_median_gatherer_sg=getSEMedian(median_gatherer_sg,1000);
        med_median_gatherer_fg=median(median_gatherer_fg);
        sem_median_gatherer_fg=getSEMedian(median_gatherer_fg,1000);
        % Sign-Rank test
        [p,h,stats] = signrank(median_gatherer_sg, median_gatherer_fg, 'alpha',0.05,'tail','right');
        disp(["p-value equals to",string(p)])
    end
    % errorbar plot
    errorbar(2,med_median_gatherer_sg,sem_median_gatherer_sg,'g','LineWidth',2,'MarkerSize',15)
    errorbar(1,med_median_gatherer_fg,sem_median_gatherer_fg,'g','LineWidth',2,'MarkerSize',15)
    %%%%%%%%%%%%%%%%%% BOXPLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % superimpose a boxplot on the violin plot
    % data=[(median_gatherer_sg)',(median_gatherer_fg)'];
    % group=({'Slow gamma','Fast gamma'});
    % boxchart(data)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    y_line = y_max + 0.03; % Add a buffer above the max value
    %%%%%%%% Horizontal plot %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if axis_ori
        % Draw the significance line
        line(x_positions, [y_line y_line], 'Color', 'k', 'LineWidth', 2);
        text(mean(x_positions), y_line + 0.005, significance, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 14, 'FontWeight', 'bold', ...
        'Rotation', -90);
        view([90 -90]); % Rotates the axes
        set(gca, 'YDir', 'normal'); % Ensure correct label orientation
    else
        % Draw the significance line
        line(x_positions, [y_line y_line], 'Color', 'k', 'LineWidth', 2);
        text(mean(x_positions), y_line + 0.005, significance, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 14, 'FontWeight', 'bold', ...
        'Rotation', 0);
    end
end