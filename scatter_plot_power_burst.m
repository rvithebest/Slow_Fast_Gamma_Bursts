function  scatter_plot_power_burst(power_gatherer_sg,power_gatherer_fg,median_gatherer_sg ...
        ,median_gatherer_fg,matched_sg_indices,matched_fg_indices,Monkey_num)
    matched_sg_power = power_gatherer_sg(matched_sg_indices);
    matched_fg_power = power_gatherer_fg(matched_fg_indices);
    matched_sg_lengths = median_gatherer_sg(matched_sg_indices);
    matched_fg_lengths = median_gatherer_fg(matched_fg_indices);
    scatter(power_gatherer_sg, median_gatherer_sg, 50, 'b', 'o', 'LineWidth', 1.5);
    hold on;
    scatter(power_gatherer_fg, median_gatherer_fg, 50, 'o', 'LineWidth', 1.5,'MarkerEdgeColor',[1 0.5 0]);
    if Monkey_num==2
      xlabel('Power difference (dB)')
    end
    ylabel('Burst lengths(s)')
    %set(gca,'FontSize',25); hold on;
    scatter(matched_sg_power, matched_sg_lengths, 50, 'bo','filled');
    scatter(matched_fg_power, matched_fg_lengths, 50, 'o', 'filled', 'MarkerFaceColor', [1 0.5 0]);
    if Monkey_num==2
        legend('Slow gamma','Fast gamma','Slow gamma(matched)','Fast gamma(matched)')
    end 
end