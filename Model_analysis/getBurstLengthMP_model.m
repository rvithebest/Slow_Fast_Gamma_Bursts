function [lengthList,freqList,timeList,gaborInfo,header,aList] = getBurstLengthMP_model(analogData,timeVals,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,gaborInfo,header,amp_temp, width_temp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get MP Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs=round(1/(timeVals(2)-timeVals(1)));

if isempty(gaborInfo)
    [B,A]=butter(4,[burstFreqRangeHz(1)-10 burstFreqRangeHz(2)+10]/(Fs/2));
    analogData1=filtfilt(B,A,analogData')';
    [gaborInfo,header] = getStochasticDictionaryMP3p1(analogData1,timeVals,maxIteration,adaptiveDictionaryParam,dictionarySize);
end
numTrials = size(gaborInfo,1);
%%%%%%%%%%%%%%%%%%%%%%%% Get Information from Atoms %%%%%%%%%%%%%%%%%%%%%%%
% Stochastic Dictionary Parameters
SCALE  =1;
FREQ   =2;
POS    =3;
MODULUS=4; %mp2tf uses modulus, although gabor uses amplitude
AMPLI  =5;
%PHASE  =6;

lengthList = cell(1,numTrials);
freqList = cell(1,numTrials);
timeList = cell(1,numTrials);
modList = cell(1,numTrials);
aList= cell(1,numTrials);
octaveToLengthMultiplier = 4/sqrt(2*pi);
for i=1:numTrials
    x = squeeze(gaborInfo(i,:,:));
    o = x(:,SCALE);
    f = x(:,FREQ);
    t = x(:,POS)+timeVals(1);
    m = x(:,MODULUS);
    a= x(:,AMPLI);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    goodAtomsInT = intersect(find(t>=stimulusPeriodS(1)),find(t<stimulusPeriodS(2))); % Find Atoms whos time center lies in the stimulus period
    goodAtomsInF = intersect(find(f>=burstFreqRangeHz(1)),find(f<burstFreqRangeHz(2))); % Find Atoms freq center lies in the band
    goodAtomsInTF = intersect(goodAtomsInF,goodAtomsInT);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mtmp = m(goodAtomsInTF);
    scaling_factor=(3.7/(mean(width_temp)^2));
    useThesePos=find(mtmp>(amp_temp(i)*scaling_factor));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(useThesePos)
        lengthList{i} = octaveToLengthMultiplier*o(useThesePos);
        freqList{i} = f(useThesePos);
        timeList{i} = t(useThesePos);
        modList{i} = m(useThesePos);
        aList{i}=a(useThesePos);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% DISPLAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayFlag
    cLims = [-3 3]; colormap jet;
    freqRange = [0 100]; timeRange = [-1 2];
    
    % Display Mean Energy and the Bands
    subplot(221);
    [meanE,freqVals]=getEnergyMP3p1(gaborInfo,header,timeVals);
    pcolor(timeVals,freqVals,log10(meanE));
    ylim(freqRange); caxis(cLims);
    shading interp;
    title('Stochastic');
    
    hold on;
    numTimeVals=length(timeVals);
    plot(timeVals,zeros(1,numTimeVals)+ burstFreqRangeHz(1),'color','k');
    plot(timeVals,zeros(1,numTimeVals)+ burstFreqRangeHz(2),'color','k');
    plot(zeros(1,length(freqVals))+ stimulusPeriodS(1),freqVals,'k--');
    plot(zeros(1,length(freqVals))+ stimulusPeriodS(2),freqVals,'k--');
    axis([timeRange freqRange]);
    
    %%%%%%%%%%%%%%%%%%%%%%% Single Trial Analysis %%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:numTrials
        disp(['Trial: ' num2str(i)]);
        
        subplot(222);
        cla;
        plot(timeVals,analogData(i,:));
        xlim(timeRange);

        subplot(224);
        cla;
        E = mp2tf(squeeze(gaborInfo(i,:,:)),header(i,:));
        pcolor(timeVals,freqVals,log10(E));
        caxis(cLims);
        shading interp;
        
        hold on;
        plot(timeVals,zeros(1,numTimeVals)+ burstFreqRangeHz(1),'color','k');
        plot(timeVals,zeros(1,numTimeVals)+ burstFreqRangeHz(2),'color','k');
        plot(zeros(1,length(freqVals))+ stimulusPeriodS(1),freqVals,'k--');
        plot(zeros(1,length(freqVals))+ stimulusPeriodS(2),freqVals,'k--');
        for k=1:length(lengthList{i})
            plot(timeList{i}(k),freqList{i}(k),'color','k','marker','o','markersize',4);
            tmpList = intersect(find(timeVals>=(timeList{i}(k)-(lengthList{i}(k))/2)),find(timeVals<= (timeList{i}(k)+(lengthList{i}(k))/2)));
            if ~isempty(tmpList)
                plot(timeVals(tmpList),zeros(1,length(tmpList))+freqList{i}(k),'color','k','linewidth',1);
            end
        end
        axis([timeRange freqRange]);
        pause;
    end
end
end
