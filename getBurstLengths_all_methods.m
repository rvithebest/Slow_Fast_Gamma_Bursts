load('LFP_synth_gamma_all_methods_elec41_ORI_45.mat');
synthetic_LFP_data=analogData_accumulator(:,1);
timeVals_all=analogData_accumulator(:,3);
thresholdFraction=0.5;
displayFlag=0;
stimulusPeriodS=[0.25 0.75];
baselinePeriodS=[-0.5 0];
burstFreqRangeHz=[20 65]; % Bursts injected in this freq range
filterOrder=4;
%%% Intializing length_gatherer variables %%%%%%%%
length_gatherer_hilbert=cell(8,1);
length_gatherer_wavelet=cell(8,1);
length_gatherer_feingold=cell(8,1);
for i=1:8
    data_temp=synthetic_LFP_data{i,1};
    timeVals=timeVals_all{i,1};
    numTrials=size(data_temp,1);
    diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz);
    thresholdFactor=(diffPower*thresholdFraction);
    %%% Hilbert Transform Method %%%%%%%
    [length_temp_hilbert,~]=getBurstLengthHilbert(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,filterOrder);
    length_temp_hilbert_all_trials=[];
    for j=1:numTrials
        length_temp_hilbert_all_trials=[length_temp_hilbert_all_trials, length_temp_hilbert{j}];
    end
    length_gatherer_hilbert{i,1}=length_temp_hilbert_all_trials;
    %%%% Wavelet Transform Method %%%%%%%
    [length_temp_wavelet,~]=getBurstLengthWavelet(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz);
    length_temp_wavelet_all_trials=[];
    for j=1:numTrials
        length_temp_wavelet_all_trials=[length_temp_wavelet_all_trials, length_temp_wavelet{j}];
    end
    length_gatherer_wavelet{i,1}=length_temp_wavelet_all_trials;
    %%%%% Feingold Method %%%%%%%
    [length_temp_feingold]=getBurstLengthFeingold(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,filterOrder);
    length_temp_feingold_all_trials=[];
    for j=1:numTrials
        length_temp_feingold_all_trials=[length_temp_feingold_all_trials, length_temp_feingold{j}];
    end
    length_gatherer_feingold{i,1}=length_temp_feingold_all_trials;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
save('LFP_synth_gamma_rest_methods_elec41_ORI_45.mat','length_gatherer_feingold','length_gatherer_wavelet','length_gatherer_hilbert');