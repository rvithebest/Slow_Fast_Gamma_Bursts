function save_MP_results(lfp,t,folder_name,indicator)
    Fs=1/(t(2)-t(1)); % Sampling frequency
    %%%%%%%%%%%% Low-pass filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Design a low-pass filter (Butterworth filter, 4th order, cutoff frequency 200 Hz)
    [b, a] = butter(4, 200/(Fs/2), 'low');
    % Apply the filter to the LFP data
    lfp_filtered = zeros(size(lfp));
    for numTrial=1:size(lfp,1)
        lfp_filtered(numTrial,:) = filtfilt(b, a, lfp(numTrial,:)')';
    end
    %%%%%%%%%%%% Decimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decimationFactor=200;
    lfp_decimated = zeros(size(lfp_filtered,1), 1024); % N=1024
    for numTrial=1:size(lfp,1)
        lfp_decimated(numTrial,:)=resample(lfp_filtered(numTrial,:)',1,decimationFactor)'; 
     end
    timeVals_decimated = downsample(t', decimationFactor)';
    %%%%%%%%%%%%%%% MP length analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    displayFlag=0;
    stimulusPeriodS=[0.25 1.25];
    baselinePeriodS=[-1 0];
    if indicator==0
        gamma_freq=[20 40]; % Slow gamma frequency
        num_iterations=50;
    else
        gamma_freq=[40 80]; % Fast gamma frequency
        num_iterations=100;
    end
    thresholdFraction=0.5;
    dict_size=2500000;
    thresholdFactor=sqrt(thresholdFraction);
    [length_temp,~,~,gabor_temp,header_temp,~]= getBurstLengthMP(lfp_decimated,timeVals_decimated,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,gamma_freq,num_iterations,0.9,dict_size,[],[]);
    length_temp_all_trials=[];
    for ii=1:length(length_temp)
         length_temp_all_trials=[length_temp_all_trials,length_temp{ii}'];
    end
    % Save results in the specified folder as MP_analysis_results.mat
    save(fullfile(folder_name, 'MP_analysis_results.mat'), 'length_temp_all_trials', 'timeVals_decimated', 'lfp_decimated', 'gabor_temp', 'header_temp', 'gamma_freq', 'num_iterations', 'dict_size','displayFlag', 'stimulusPeriodS', 'baselinePeriodS', 'indicator','thresholdFactor');
    disp('MP analysis results saved successfully.');
end