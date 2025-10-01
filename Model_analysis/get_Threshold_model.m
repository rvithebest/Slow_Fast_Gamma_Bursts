function [Thresh_vals, Width_vals] = get_Threshold_model(analogData, timeVals, timePeriod, indicator)
    if indicator==1
        gammaRange = [20 40]; % Slow gamma range- proxy for baseline
    elseif indicator==2
        gammaRange = [40 80]; % Fast gamma range- proxy for baseline
    end
    % Calculate indices dynamically
    fs = 1 / (timeVals(2) - timeVals(1)); % Sampling frequency
    %timePeriod of interest
    timeIndices = find(timeVals >timePeriod(1) & timeVals < timePeriod(2));    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make power-spectral density using multitaper chronux package
    TW = 2; % time-bandwidth product
    K = 3; % number of tapers
    params.Fs = fs;
    params.tapers = [TW K];
    params.fpass = [0 250];
    params.trialave = 0;
    params.pad = -1; % Disable zero-padding
    params.err = [2 0.05];
    [S, f] = mtspectrumc(analogData(:, timeIndices)', params);
    S=S';
    % Select the frequency range of interest - gammaRange
    f_gamma = f(f >= gammaRange(1) & f <= gammaRange(2));
    S_gamma = S(:,f >= gammaRange(1) & f<= gammaRange(2));
    Peak_power = log10(max(S_gamma,[], 2));
    Norm_factor=(log10(S_gamma(:,1))+log10(S_gamma(:,end)))/2;
    % Thresh_vals=10*(Peak_power-Norm_factor);
    % Thresh_vals=sqrt(max(S_gamma,[], 2));
    scale_factor=10;
    Thresh_vals=Peak_power+scale_factor;
    Half_power=((Peak_power-Norm_factor)/2)+Norm_factor;
    Width_vals=nan(size(Peak_power));
    S_gamma=log10(S_gamma);
    for i=1:size(S_gamma,1)
        freq_above_half=find(S_gamma(i,:)>=Half_power(i));
        if ~isempty(freq_above_half)
            Width_vals(i)=f_gamma(freq_above_half(end))-f_gamma(freq_above_half(1));
        else
            Width_vals(i)=0;
        end
    end
end
