function Power = get_Power(analogData, timeVals, timePeriod, gammaRange)    
    % Calculate indices dynamically
    fs = 1 / (timeVals(2) - timeVals(1)); % Sampling frequency
    %timePeriod of interest
    timeIndices = find(timeVals >timePeriod(1) & timeVals < timePeriod(2));    
    %%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    % Exclude the 50 Hz point from the gamma range using a tolerance
    % tol = 1e-6; % Tolerance for floating-point comparison
    % exclude_50Hz_index = find(abs(f_gamma - 50) < tol);
    % if gammaRange(1)>44
    %     if (exclude_50Hz_index==[])
    %         error('Stop')
    %     end
    % end
    % if ~isempty(exclude_50Hz_index)
    %     f_gamma(exclude_50Hz_index) = [];
    %     S_gamma = S(f >= gammaRange(1) & f <= gammaRange(2), :);
    %     S_gamma(exclude_50Hz_index, :) = [];
    % else
    %     S_gamma = S(f >= gammaRange(1) & f<= gammaRange(2), :);
    % end
    S_gamma = S(f >= gammaRange(1) & f<= gammaRange(2), :);
    % TAKE THE MEAN OF THE S_gamma in 1st dimension (across frequency range)
    %Power=S_gamma;
    Power = mean(S_gamma, 1);
    % now , take average of the S_gamma in the 2nd dimension(across trials)
    % Power=mean(S_gamma, 2);
end
