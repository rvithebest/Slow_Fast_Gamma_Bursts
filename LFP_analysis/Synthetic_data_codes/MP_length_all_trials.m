length_gatherer = [];
length_gatherer_beta=[];
length_gatherer_gamma=[];
length_gatherer_delta=[];
length_gatherer_MP=[];
length_gatherer_MP_beta=[];
length_gatherer_MP_gamma=[];
length_gatherer_MP_delta=[];
length_gatherer_hilbert=[];
length_gatherer_hilbert_beta=[];
length_gatherer_hilbert_gamma=[];
length_gatherer_hilbert_delta=[];
length_gatherer_wavelet=[];
length_gatherer_wavelet_beta=[];
length_gatherer_wavelet_gamma=[];
length_gatherer_wavelet_delta=[];
length_gatherer_feingold=[];
length_gatherer_feingold_beta=[];
length_gatherer_feingold_gamma=[];
length_gatherer_feingold_delta=[];
thresholdFraction=0.5;
displayFlag=0;
stimulusPeriodS=[0.25 0.75];
baselinePeriodS=[-0.5 0];
%%% narrow range is taken for OMP-GEAR(separately for slow and fast
%%% ranges if required
burstFreqRangeHz=[20 65];
numTrials=size(analogData,1);
diffPower = getChangeInPower(analogData,timeVals,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz);
% thresholdFactor=0.38*sqrt(diffPower);-optimized threshold for
% synthetic_data_fourier(Decimated EEG data-good performance
thresholdFactor=sqrt(diffPower*thresholdFraction);
% [length_measured,~,~,~] = getBurstLength_all(analogData(i,:),timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,50,0.9,[],[],[],"OMP-MAGE");
[length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,1500000,[],[],'OMP-GEAR');
% length_gatherer = [length_gatherer, length_measured{1,1}', length_measured{1,2}', length_measured{1,3}'];   
for ii=1:numTrials
     length_gatherer = [length_gatherer, length_measured{1,ii}'];
end
thresholdFactor=sqrt(diffPower*0.25);
[length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,1500000,gaborInfo,header,'OMP-GEAR');
for ii=1:numTrials
    length_gatherer_beta = [length_gatherer_beta, length_measured{1,ii}'];
end
thresholdFactor=sqrt(diffPower*0.75);
[length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,1500000,gaborInfo,header,'OMP-GEAR');
for ii=1:numTrials
    length_gatherer_gamma = [length_gatherer_gamma, length_measured{1,ii}'];
end
thresholdFactor=sqrt(diffPower*0.14);
[length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,1500000,gaborInfo,header,'OMP-GEAR');
for ii=1:numTrials
    length_gatherer_delta = [length_gatherer_delta, length_measured{1,ii}'];
end
gaborInfo_current=gaborInfo;
% MP
thresholdFactor=sqrt(diffPower*thresholdFraction);
[length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,2500000,[],[],'MP');
for ii=1:numTrials
    length_gatherer_MP = [length_gatherer_MP, length_measured{1,ii}'];
end
thresholdFactor=sqrt(diffPower*0.25);
[length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,2500000,gaborInfo,header,'MP');
for ii=1:numTrials
    length_gatherer_MP_beta = [length_gatherer_MP_beta, length_measured{1,ii}'];
end
thresholdFactor=sqrt(diffPower*0.75);
[length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,2500000,gaborInfo,header,'MP');
for ii=1:numTrials
    length_gatherer_MP_gamma = [length_gatherer_MP_gamma, length_measured{1,ii}'];
end
thresholdFactor=sqrt(diffPower*0.14);
[length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,2500000,gaborInfo,header,'MP');
for ii=1:numTrials
    length_gatherer_MP_delta = [length_gatherer_MP_delta, length_measured{1,ii}'];
end
gaborInfo_MP=gaborInfo;
%%%% Hilbert transform %%%%%%%%%%%%%%%
thresholdFactor=(diffPower*thresholdFraction);
length_measured = getBurstLengthHilbert(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,4);
for ii=1:numTrials
    length_gatherer_hilbert = [length_gatherer_hilbert, length_measured{1,ii}];
end
thresholdFactor=(diffPower*0.25);
length_measured = getBurstLengthHilbert(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,4);
for ii=1:numTrials
    length_gatherer_hilbert_beta = [length_gatherer_hilbert_beta, length_measured{1,ii}];
end
thresholdFactor=(diffPower*0.75);
length_measured = getBurstLengthHilbert(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,4);
for ii=1:numTrials
    length_gatherer_hilbert_gamma = [length_gatherer_hilbert_gamma, length_measured{1,ii}];
end
thresholdFactor=(diffPower*0.14);
length_measured = getBurstLengthHilbert(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,4);
for ii=1:numTrials
    length_gatherer_hilbert_delta = [length_gatherer_hilbert_delta, length_measured{1,ii}];
end
%%%% Wavelet transform %%%%%%%%%%%%%%%%%%%%%%
thresholdFactor=(diffPower*thresholdFraction);
[length_temp_wavelet,~]= getBurstLengthWavelet(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz);
for ii=1:numTrials
    length_gatherer_wavelet = [length_gatherer_wavelet, length_temp_wavelet{1,ii}];
end
thresholdFactor=(diffPower*0.25);
[length_temp_wavelet,~]= getBurstLengthWavelet(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz);
for ii=1:numTrials
    length_gatherer_wavelet_beta = [length_gatherer_wavelet_beta, length_temp_wavelet{1,ii}];
end
thresholdFactor=(diffPower*0.75);
[length_temp_wavelet,~]= getBurstLengthWavelet(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz);
for ii=1:numTrials
    length_gatherer_wavelet_gamma = [length_gatherer_wavelet_gamma, length_temp_wavelet{1,ii}];
end
thresholdFactor=(diffPower*0.14);
[length_temp_wavelet,~]= getBurstLengthWavelet(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz);
for ii=1:numTrials
    length_gatherer_wavelet_delta = [length_gatherer_wavelet_delta, length_temp_wavelet{1,ii}];
end
%%%% Feingold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresholdFactor=(diffPower*thresholdFraction);
[length_temp_feingold]= getBurstLengthFeingold(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,4);
for ii=1:numTrials
    length_gatherer_feingold = [length_gatherer_feingold, length_temp_feingold{1,ii}];
end
thresholdFactor=(diffPower*0.25);
[length_temp_feingold]= getBurstLengthFeingold(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,4);
for ii=1:numTrials
    length_gatherer_feingold_beta = [length_gatherer_feingold_beta, length_temp_feingold{1,ii}];
end
thresholdFactor=(diffPower*0.75);
[length_temp_feingold]= getBurstLengthFeingold(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,4);
for ii=1:numTrials
    length_gatherer_feingold_gamma = [length_gatherer_feingold_gamma, length_temp_feingold{1,ii}];
end
thresholdFactor=(diffPower*0.14);
[length_temp_feingold]= getBurstLengthFeingold(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,4);
for ii=1:numTrials
    length_gatherer_feingold_delta = [length_gatherer_feingold_delta, length_temp_feingold{1,ii}];
end
