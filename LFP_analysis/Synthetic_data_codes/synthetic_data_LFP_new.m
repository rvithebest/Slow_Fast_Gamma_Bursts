function [analogData,timeVals,analogData0] = synthetic_data_LFP(burstLen)
    if nargin<1
        burstLen=0.4;
    end
    %default parameters for the stimulation
    cvAmp = 0.1;                
    displayFlag = 1;                         
    stimulusPeriod = [0.25 0.75];     
    %%%%%% BROAD GAMMA RANGE %%%%%%%%%%%%%%%%%%%%
    gammaRange = [20 65];          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % numBurstsPerTrial=1;           
    % if burstLen<0.15
    numBurstsPerTrial=[];
    if burstLen>0.25
        numBurstsPerTrial=1; %makes sure that burst is injected
    end
    % end
    % Get Real Data
    % We select electrode number 8
    % ORI=157.5 deg and SF=1 cpd and 4 cpd 
    % Corresponding trials are taken
    load("Decimated_8_LFP_data_alpa_H\deci_elec_41.mat")
    load("Decimated_8_LFP_data_alpa_H\timeVals_decimated.mat")
    timeVals=timeVals_decimated;
    % load the parameterCombinations.mat file which gives the 'SF' and 'ORI' information
    load("alpaH_info\parameterCombinations.mat");
    %badtrials file is loaded
    load("alpaH_info\badTrials.mat");
    trial_pos=(parameterCombinations{:,:,:,2,9});
    % trial_pos_2=(parameterCombinations{:,:,:,3,3});
    % trial_pos_3=(parameterCombinations{:,:,:,4,3});
    % slow_gamma_pos=[trial_pos_1,trial_pos_2,trial_pos_3];
    trials_selected=setdiff(trial_pos,badTrials);
    LFP_data=analogDataDecimated(trials_selected,:);
    bl= LFP_data(:,84:211);
    st = LFP_data(:,274:401);
    fft_bl=fft(bl,[],2);
    ang_fft_bl=angle(fft_bl);
    abs_fft_bl=abs(fft_bl);
    avg_fft_bl=mean(abs_fft_bl);
    freq_axis_original = linspace(0, 1-(1/128), 128)*250;
    % Target points (x-coordinates) for interpolation
    freq_axis_interpolated = linspace(0, 1-(1/512), 512)*250;
    % Preallocate interpolated data array
    log_fft_interpolated_magnitude = zeros(1, length(freq_axis_interpolated));
     % Interpolate the data
    log_fft_bl = log10(avg_fft_bl(:, 1:65));
    log_fft_interpolated_magnitude(:, 1) = log_fft_bl(:, 1);
    % reduces the effct of DC component- emperically obtained
    log_fft_interpolated_magnitude(:, 2:257) = interp1(freq_axis_original(1:65), log_fft_bl, freq_axis_interpolated(2:257), 'linear');
    log_fft_interpolated_magnitude(:, 258:512) = log_fft_interpolated_magnitude(:, 256:-1:2);
    fft_interpolated_magnitude = 10.^(log_fft_interpolated_magnitude);
    %normalize the fft as the number of points changes
    normalization_factor=(sqrt(512/128));
    fft_interpolated_magnitude = fft_interpolated_magnitude * normalization_factor;
    % Generate random phases and ensure Hermitian symmetry
    num_trials = size(LFP_data,1);
    interpolated_points = size(fft_interpolated_magnitude, 2);
    half_points = (interpolated_points / 2) + 1;
    random_phase = zeros(num_trials, interpolated_points);
    random_phase(:, 1) = 0;
    random_phase(:, half_points) = 0;
    random_phase(:, 2:half_points-1) = 2 * pi * rand(num_trials, half_points - 2) - pi;
    
    % keep  the phase of the original signal for the original points-(5,9,13,....,509)
    for i = 1:num_trials
        k=2;
        for j=5:4:509
            if (k==65)
                break
            end
            random_phase(i, j) = ang_fft_bl(i, k)-pi;
            k=k+1;
        end
    end
    random_phase(:, half_points + 1:end) = -random_phase(:, half_points - 1:-1:2);
    % Generate the complex FFT values
    fft_interpolated_magnitude = repmat(fft_interpolated_magnitude, num_trials, 1);
    fft_interpolated_complex = fft_interpolated_magnitude .* exp(1i * random_phase);
    
    spontaneous_signal = ifft(fft_interpolated_complex, [], 2);
    
    Fs = round(1/(timeVals(2)-timeVals(1)));
    freqVals=freq_axis_original;
    stPos = 274:401;
    blPos=84:211;
    fPos = intersect(find(freqVals>=gammaRange(1)),find(freqVals<gammaRange(2)));
    numGammaPoints = length(fPos);
    numTrials=size(LFP_data,1);
    % Get Mean FFTs
    mFFTst = mean(abs(fft(st(1:numTrials,:),[],2)));
    mFFTbl = mean(abs(fft(bl(1:numTrials,:),[],2)));
    
    % Get parameters for simulations
    meanGammaAmps = (mFFTst(fPos)- mFFTbl(fPos));
    % burstLen=0.3;
    burstSignal = zeros(numTrials,length(timeVals));
    for i=1:numTrials
        burstCenterTimeList = getBurstTimes(numBurstsPerTrial,stimulusPeriod,burstLen);
        numBursts = length(burstCenterTimeList);
        if numBursts>0
            for j=1:numBursts
                burstCenterTime = burstCenterTimeList(j);
                fIndex = randi(numGammaPoints); % uniformly within gamma
                burstCenterFreq = freqVals(fPos(fIndex));
                burstAmplitude = meanGammaAmps(fIndex)*(1+cvAmp*randn);
                burst = burstAmplitude*cos(2*pi*burstCenterFreq*timeVals+2*pi*rand).*exp(-((timeVals-burstCenterTime)/(sqrt(2)*(burstLen/4))).^2);
                burstSignal(i,:) = burstSignal(i,:)+burst;
            end
        end
    end
    fft_spontaneous=fft(spontaneous_signal,[],2);
    fft_LFP=fft(LFP_data,[],2);
    fft_spontaneous(:,1)=fft_LFP(:,1);% matches the DC
    spontaneous_signal=ifft(fft_spontaneous,[],2);
    
    analogData0 = spontaneous_signal;
    analogData = analogData0 + burstSignal;
    fft_synth_st=mean(abs(fft(analogData(:, stPos), [], 2)));
    fft_synth_bl=mean(abs(fft(analogData(:, blPos), [], 2)));
    fft_spon_st= mean(abs(fft(analogData0(:,stPos), [], 2)));
    fft_spon_bl= mean(abs(fft(analogData0(:,blPos), [], 2)));
    % Rescale the burst sizes such that gamma band has the same energy
    gammaPowerST = sum(mFFTst(fPos).^2);
    gammaPowerSynth = sum(fft_synth_st(fPos).^2);
    scalingFactor = (0.7)*sqrt(gammaPowerST/gammaPowerSynth);
    % scalingFactor=scalingFactor*0.5;
    analogData = analogData0 + (scalingFactor)*burstSignal;
    if displayFlag
        figure;
        plot(freq_axis_original,log10(avg_fft_bl), 'ro', 'DisplayName', 'Original 128 points-baseline');
        hold on;
        plot(freq_axis_interpolated, log_fft_interpolated_magnitude(1, :), 'b*', 'DisplayName', 'Interpolated 512 points (linear)');
        legend;
        xlabel('Frequency');
        ylabel('Log(Amplitude)');
        title('FFT Interpolation from 128 to 512 Points (Linear, Trial avergaed)');
        synthColorName = [285, 192, 203] / 285;
        figure;
        % plot(freqVals, 10*log10(mFFTst./mFFTbl), 'k'); hold on;
        params.Fs=250;
        params.fpass=[0 100];
        params.tapers=[1 1];
        params.trialave=1;
        params.pad=-1;
        params.error=[2 0.05];
        [S_stim_o,f_stim_o]=mtspectrumc(st',params);
        [S_bl_o,f_bl_o]=mtspectrumc(bl',params);
        plot(f_stim_o,10*log10(S_stim_o./S_bl_o),'LineWidth',2,'Color','r');
        hold on;
        [S,f]=mtspectrumc(analogData(:, stPos)',params);
        [S0,f0]=mtspectrumc(analogData0(:,stPos)',params);
        [S_bl,f_bl]=mtspectrumc(analogData(:, blPos)',params);
        [S0_bl,f0_bl]=mtspectrumc(analogData0(:,blPos)',params);
        plot(f,10*log10(S./S_bl),'LineWidth',2,'Color','k');
        hold on;
        plot(f0,10*log10(S0./S0_bl),'LineWidth',2,'Color','g');
        ylabel('Power(dB)');
        % set(gca,'FontSize',18);
        xlabel('Frequency (Hz)')
        % plot(freqVals, log10(mFFTbl), 'g');
        % plot(freq_axis_interpolated, log10(mean(abs(fft(LFP_data(:,:), [], 2)))), 'color', 'c');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fft_synth_st=mean(abs(fft(analogData(:, stPos), [], 2)));
        % fft_synth_bl=mean(abs(fft(analogData(:, blPos), [], 2)));
        % fft_spon_st= mean(abs(fft(analogData0(:,stPos), [], 2)));
        % fft_spon_bl= mean(abs(fft(analogData0(:,blPos), [], 2)));
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot(freqVals, 10*log10(fft_synth_st./fft_synth_bl), 'color', 'r');
        % plot(freqVals, 10*log10(fft_spon_st./fft_spon_bl), 'color', 'b');
        xlim([0 100]);
        % plot(freq_axis_interpolated, log10(mean(abs(fft(analogData0(:,:), [], 2)))), 'color', 'r');
        % calculate the energy in the signal of each type (average over trials) in frequency domain
        energySynth=mean(sum(abs(fft(analogData(:,stPos), [], 2)).^2,2));
        energySynth0=mean(sum(abs(fft(analogData0(:,stPos), [], 2)).^2,2));
        energy_spontaneous=mean(sum(abs(fft(analogData0(:,:), [], 2)).^2,2));
        energyEEG=mean(sum(abs(fft(LFP_data(:,:), [], 2)).^2,2));
        % title(['Synthetic Signal Energy:',num2str(energySynth),'; Spontaneous Signal Energy:',num2str(energySynth0),'; EEG Energy:',num2str(energyEEG),'; Spontaneous Energy_comp:',num2str(energy_spontaneous)]);
        % legend('stimulus', 'baseline','injected-st','spontaneous-st');
        legend("Original","Synthetic","Spontaneous");
        figure;
        plot(timeVals,analogData(1,:)); hold on;
        plot(timeVals,analogData0(1,:)); hold on;
        plot(timeVals,LFP_data(1,:));
        %calculate the energy in the signal of each type (average over trials)
        energySynth=mean(sum(analogData(:,:).^2,2));
        energySynth0=mean(sum(analogData0(:,:).^2,2));
        energyEEG=mean(sum(LFP_data(:,:).^2,2));
        title(['Synthetic Signal Energy:',num2str(energySynth),'; Spontaneous Signal Energy:',num2str(energySynth0),'; EEG Energy:',num2str(energyEEG)]);
        xlabel('Time(seconds)')
        legend('Signal-after injection', 'sponataneous', 'original EEG');
    end

end
function burstCenterTimeList = getBurstTimes(numBurstsPerTrial,stimulusPeriod,burstLen)

if isempty(numBurstsPerTrial)
    numBursts0 = poissrnd(diff(stimulusPeriod)/burstLen); % Intial number of bursts
else
    numBursts0 = numBurstsPerTrial;
end

burstCenterTimeRange=[(stimulusPeriod(1)+(burstLen/2)) (stimulusPeriod(2)-(burstLen/2))];    % Limits on the Mean of the Gabor signal to be injected

if numBursts0==0
    burstCenterTimeList = [];
elseif numBursts0==1
    burstCenterTimeList = burstCenterTimeRange(1) + diff(burstCenterTimeRange)*rand; %  take one time uniformly within allowed time range
else
    burstCenterTimeList0 = burstCenterTimeRange(1) + diff(burstCenterTimeRange)*rand(1,numBursts0); % create a list of numBurst0 burst times
    sortedBurstCenterTimeList0 = sort(burstCenterTimeList0); % Sort this list
    
    tPointer = sortedBurstCenterTimeList0(1);
    burstCenterTimeList(1) = tPointer; % first burst time is the earliest one
    
    while(tPointer<burstCenterTimeRange(2))
        nextGoodEntryPos = find(sortedBurstCenterTimeList0>tPointer+burstLen,1);
        if isempty(nextGoodEntryPos)
            tPointer=burstCenterTimeRange(2);
        else
            tPointer = sortedBurstCenterTimeList0(nextGoodEntryPos);
            burstCenterTimeList = cat(2,burstCenterTimeList,tPointer);
        end
    end
end
end