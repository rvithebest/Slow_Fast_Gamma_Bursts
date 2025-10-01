function dEndB=normalize_TF(raw_power_spetrogram,timeVals,baselineS)
    % First we calculate the log of the raw power spectrogram given as input
    threshold=10^(-16);
    % we take a threshold so that calculation of log does not give error
    raw_power_spetrogram(raw_power_spetrogram<threshold)=threshold;
    logMeanEn = log10(raw_power_spetrogram); % log of raw power spectrogram
    % Change in power from baseline
    blL = find(timeVals>=baselineS(1),1);           % lower index of baseline
    blU = find(timeVals<baselineS(2),1,'last');     % upper index of baseline
    baselineEn=mean(logMeanEn(blL:blU,:),1);        % baseline TF Energy Matrix
    dEndB = 10*(logMeanEn'-repmat(baselineEn,length(timeVals),1)'); % change in dB
    % We multipy by 10 so as to get the change in dB
end
    