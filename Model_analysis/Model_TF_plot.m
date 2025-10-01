function Model_TF_plot(LFP_signal,timeVals,idx,plotHandles,freq_range)
    Fs=1/(timeVals(2)-timeVals(1));
    stimulusPeriodS=[0.25 1.25];
    baselinePeriodS=[-1 0];
    params.Fs=Fs;
    params.fpass=[0 101];
    params.tapers=[1 1];
    params.trialave=1;
    params.pad=-1;
    params.err=[2 0.05];
    movingWin=[0.25 0.025];
    [S_temp,t_temp,f_temp]= mtspecgramc(LFP_signal',movingWin,params);
    t_temp = t_temp+timeVals(1)-(1/Fs); % Center the times with respect to the stimulus onset time
    % [S_dB]=normalize_TF(S_temp,t_temp,baselinePeriodS);
    S_dB=log10(S_temp)';
    subplot(plotHandles(idx,1));
    hold on;
    pcolor(t_temp,f_temp,S_dB);
    shading interp;
    colormap('jet');
    x=min(S_dB(:));
    y=max(S_dB(:));
    clim([-8 -2]);
    % draw a horizontal lines denoting the frequency range of interest
    % yline([freq_range(1) freq_range(1)],'Color','k','LineStyle','--','LineWidth',2);
    % yline([freq_range(2) freq_range(2)],'Color','k','LineStyle','--','LineWidth',2);
    c=colorbar;
    c.Position = c.Position + [-0.15 0 0.003 0];
    c.Ticks  = [-8 -5 -2];
    c.FontSize = 18;
    xlim([-0.5 1.8]);
    % xlim([-0.5 1.6]);
    ylim([0 100]);
    if idx==2
        xlabel("Time (s)")
    end
    ylabel("Frequency (Hz)")
    hold on;
    % draw a black lines representing the stimulus period (vertical dashed lines)
    line([stimulusPeriodS(1) stimulusPeriodS(1)], ylim, 'Color', 'k', 'LineStyle', '--','LineWidth',2);
    line([stimulusPeriodS(2) stimulusPeriodS(2)], ylim, 'Color', 'k', 'LineStyle', '--','LineWidth',2);
end