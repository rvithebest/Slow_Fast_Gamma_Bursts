function plot_burst_lengths(length_measured_sg,time_measured_sg,trial,color_plot)
    burst_length_sg=length_measured_sg{1,trial};
    burst_time_sg=time_measured_sg{1,trial};
    % draw the bursts, with burst (burst_length_sg(ii)) centered at time=burst_time_sg(ii)
    hold on; % Keep all bursts on the same plot
    for ii = 1:length(burst_length_sg)
        if burst_length_sg(ii)>0.8
            continue
        end
        % Calculate start and end points of the burst
        start_time = burst_time_sg(ii) - burst_length_sg(ii)/2;
        end_time = burst_time_sg(ii) + burst_length_sg(ii)/2;
        % Plot the burst as a horizontal line at y = 0
        plot([start_time, end_time], [0, 0], 'LineWidth', 2,'Color','k'); % Line at y=0
        % plot burst center as a red dot
        plot(burst_time_sg(ii), 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', color_plot);
    end
    % draw two dotted lines to indicate the start and end of the stimulus period - [0.25 0.75]- color magenta
    line([0.25 0.25], [-170 160], 'Color', 'm', 'LineStyle', '--','LineWidth',2);
    line([0.75 0.75], [-170 160], 'Color', 'm', 'LineStyle', '--','LineWidth',2);
    % baseline period is [-0.5 0]- brown color
    line([-0.5 -0.5], [-100 100], 'Color', [0.23,0.51,0], 'LineStyle', '--','LineWidth',2);
    line([0 0], [-100 100], 'Color', [0.23,0.51,0], 'LineStyle', '--','LineWidth',2);
    xlim([0 1])
    % xlabel('Time (s)');
    %ylabel('Amplitude (uV)');
end