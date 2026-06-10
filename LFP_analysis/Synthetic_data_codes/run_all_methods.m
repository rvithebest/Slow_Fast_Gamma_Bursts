function [length_gatherer,length_gatherer_tf_l,length_gatherer_tf_h,...
    length_gatherer_MP,length_gatherer_MP_tf_l,length_gatherer_MP_tf_h,...
    length_gatherer_hilbert,length_gatherer_hilbert_tf_l,length_gatherer_hilbert_tf_h,...
    length_gatherer_wavelet,length_gatherer_wavelet_tf_l,length_gatherer_wavelet_tf_h] = run_all_methods(analogData,timeVals)
    % Wavelet transform not used in the manuscript (length biased by frequency), but code is here for completeness.
    length_gatherer = [];
    length_gatherer_tf_l= [];
    length_gatherer_tf_h= []
    length_gatherer_MP = [];
    length_gatherer_MP_tf_l = [];
    length_gatherer_MP_tf_h = [];
    length_gatherer_hilbert=[];
    length_gatherer_hilbert_tf_l=[];
    length_gatherer_hilbert_tf_h=[];
    length_gatherer_wavelet=[];
    length_gatherer_wavelet_tf_l=[];
    length_gatherer_wavelet_tf_h=[];
    thresholdFraction=0.5;
    tf_l=0.25;
    tf_h=0.75;
    displayFlag=0;
    stimulusPeriodS=[0.25 0.75];
    baselinePeriodS=[-0.5 0];
    burstFreqRangeHz=[20 65];
    numTrials=size(analogData,1);
    diffPower = getChangeInPower(analogData,timeVals,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz);
    thresholdFactor=sqrt(diffPower*thresholdFraction); % default threhsold
    [length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,1500000,[],[],'OMP-GEAR'); 
    for ii=1:numTrials
         length_gatherer = [length_gatherer, length_measured{1,ii}'];
    end
    thresholdFactor=sqrt(diffPower*tf_l);
    [length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,1500000,gaborInfo,header,'OMP-GEAR');
    for ii=1:numTrials
        length_gatherer_tf_l = [length_gatherer_tf_l, length_measured{1,ii}'];
    end
    thresholdFactor=sqrt(diffPower*tf_h);
    [length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,1500000,gaborInfo,header,'OMP-GEAR');
    for ii=1:numTrials
        length_gatherer_tf_h = [length_gatherer_tf_h, length_measured{1,ii}'];
    end
    gaborInfo_current=gaborInfo;
    % MP
    thresholdFactor=sqrt(diffPower*thresholdFraction);
    [length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,2500000,[],[],'MP');
    for ii=1:numTrials
        length_gatherer_MP = [length_gatherer_MP, length_measured{1,ii}'];
    end
    thresholdFactor=sqrt(diffPower*tf_l);
    [length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,2500000,gaborInfo,header,'MP');
    for ii=1:numTrials
        length_gatherer_MP_tf_l = [length_gatherer_MP_tf_l, length_measured{1,ii}'];
    end
    thresholdFactor=sqrt(diffPower*tf_h);
    [length_measured,freqList,timeList,gaborInfo,header,modList]= getBurstLength_all(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,120,0.9,2500000,gaborInfo,header,'MP');
    for ii=1:numTrials
        length_gatherer_MP_tf_h = [length_gatherer_MP_tf_h, length_measured{1,ii}'];
    end
    gaborInfo_MP=gaborInfo;
    %%%%%%% Hilbert transform %%%%%%%%%%%%%%%
    thresholdFactor=(diffPower*thresholdFraction);
    length_measured = getBurstLengthHilbert(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,4);
    for ii=1:numTrials
        length_gatherer_hilbert = [length_gatherer_hilbert, length_measured{1,ii}];
    end
    thresholdFactor=(diffPower*tf_l);
    length_measured = getBurstLengthHilbert(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,4);
    for ii=1:numTrials
        length_gatherer_hilbert_tf_l = [length_gatherer_hilbert_tf_l, length_measured{1,ii}];
    end
    thresholdFactor=(diffPower*tf_h);
    length_measured = getBurstLengthHilbert(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,4);
    for ii=1:numTrials
        length_gatherer_hilbert_tf_h = [length_gatherer_hilbert_tf_h, length_measured{1,ii}];
    end
    %%%%%%%%%% Wavelet transform %%%%%%%%%%%%%%%%%%%%%%
    thresholdFactor=(diffPower*thresholdFraction);
    [length_temp_wavelet,~]= getBurstLengthWavelet(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz);
    for ii=1:numTrials
        length_gatherer_wavelet = [length_gatherer_wavelet, length_temp_wavelet{1,ii}];
    end
    thresholdFactor=(diffPower*tf_l);
    [length_temp_wavelet,~]= getBurstLengthWavelet(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz);
    for ii=1:numTrials
        length_gatherer_wavelet_tf_l = [length_gatherer_wavelet_tf_l, length_temp_wavelet{1,ii}];
    end
    [length_temp_wavelet,~]= getBurstLengthWavelet(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz);
    for ii=1:numTrials
        length_gatherer_wavelet_tf_h = [length_gatherer_wavelet_tf_h, length_temp_wavelet{1,ii}];
    end
end