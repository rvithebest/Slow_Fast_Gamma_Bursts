function [Norm_Power] = get_Power_model(analogData, timeVals, timePeriod, indicator)
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
    params.trialave = 1;
    params.pad = -1; % Disable zero-padding
    params.err = [2 0.05];
    [S, f] = mtspectrumc(analogData(:, timeIndices)', params);
    % Select the frequency range of interest - gammaRange
    f_gamma = f(f >= gammaRange(1) & f <= gammaRange(2));
    S_gamma = S(f >= gammaRange(1) & f<= gammaRange(2), :);
    % TAKE THE MEAN OF THE S_gamma in 1st dimension (across frequency range)
    % Power = log10(mean(S_gamma,1));
    % Power=S_gamma;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Peak_Amp=sqrt(max(S_gamma));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Peak_power = log10(max(S_gamma));
    % On either side of the peak, we take the mean of the first and last value (first non-decreasing across f(entire range))
    % Iterate from the idx of the peak to the left and right to find the first non-decreasing value
    % Find the idx within S (complete frequency range)
    % [~, peak_idx_gamma] = max(S_gamma); % this gives index within gamma band only
    % peak_idx = peak_idx_gamma + find(f >= gammaRange(1), 1) - 1; % if you need full-S indexing
    % left_idx = peak_idx;
    % right_idx = peak_idx;
    % S=log10(S);
    % Find the first non-decreasing value to the left
    % while left_idx > 1 && S(left_idx - 1) <= S(left_idx)
    %     left_idx = left_idx - 1;
    % end
    % Find the first non-decreasing value to the right
    % while right_idx < length(S) && S(right_idx + 1) <= S(right_idx)
    %     right_idx = right_idx + 1;
    % end
    % disp(['Peak power at frequency ' num2str(f(peak_idx)) ' Hz: ' num2str(Peak_power)]);
    % disp(['Left index: ' num2str(f(left_idx)) ', Right index: ' num2str(f(right_idx))]);
    % % Now we have the left and right indices, we can calculate the mean
    % % Norm_factor = ((S(left_idx) + S(right_idx)) / 2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Norm_factor=(log10(S_gamma(1))+log10(S_gamma(end)))/2;
    Norm_Power=10*(Peak_power-Norm_factor);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Norm_Power=Peak_power;
end
