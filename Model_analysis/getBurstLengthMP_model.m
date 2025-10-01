function [lengthList,freqList,timeList,gaborInfo,header,aList] = getBurstLengthMP_model(analogData,timeVals,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,gaborInfo,header,amp_temp, width_temp)

if ~exist('displayFlag','var');         displayFlag=1;                  end
if ~exist('stimulusPeriodS','var');     stimulusPeriodS=[0.5 1.5];      end
if ~exist('baselinePeriodS','var');     baselinePeriodS=[-1 0];         end
if ~exist('burstFreqRangeHz','var');    burstFreqRangeHz=[40 60];       end
if ~exist('maxIteration','var');        maxIteration=50;                end
if ~exist('adaptiveDictionaryParam','var');adaptiveDictionaryParam=0.9; end
if ~exist('dictionarySize','var');      dictionarySize=[];              end
if ~exist('gaborInfo','var');           gaborInfo=[];                   end

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
mBL=getMeanBaseline(gaborInfo,timeVals(1),baselinePeriodS,burstFreqRangeHz);
for i=1:numTrials
    x = squeeze(gaborInfo(i,:,:));
    o = x(:,SCALE);
    f = x(:,FREQ);
    t = x(:,POS)+timeVals(1);
    m = x(:,MODULUS);
    a= x(:,AMPLI);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reject spurious burst detection here itself
    % length_list_temp= octaveToLengthMultiplier*o;
    % reject_idx= find(length_list_temp>2);
    % o(reject_idx)=[];
    % f(reject_idx)=[];
    % t(reject_idx)=[];
    % m(reject_idx)=[];
    % a(reject_idx)=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    goodAtomsInT = intersect(find(t>=stimulusPeriodS(1)),find(t<stimulusPeriodS(2))); % Find Atoms whos time center lies in the stimulus period
    goodAtomsInF = intersect(find(f>=burstFreqRangeHz(1)),find(f<burstFreqRangeHz(2))); % Find Atoms freq center lies in the band
    goodAtomsInTF = intersect(goodAtomsInF,goodAtomsInT);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mtmp = m(goodAtomsInTF);
    % a_select=a(goodAtomsInTF);
    % scaling_factor=28;
    % 40*0.75=30
    % scaling_factor=(30/(width_temp(i)*(burstFreqRangeHz(2)-burstFreqRangeHz(1))));% input-8: 0.15 (1,1); 0.01-0.03 (for higher theta)
    % scaling_factor=(0.6/mean(width_temp));% input-8: 0.15 (1,1); 0.01-0.03 (for higher theta)
    scaling_factor=(3.7/(mean(width_temp)^2));% input-8: 0.15 (1,1); 0.01-0.03 (for higher theta)
    useThesePos=find(mtmp>(amp_temp(i)*scaling_factor));
    % s_fact=0.1;
    % useThesePos=find(mtmp>(min(mtmp)*((10^(amp_temp/10))*s_fact)));
    % useThesePos=find(mtmp>(mBL*(sqrt(10^(amp_temp/10))*s_fact)));
    % useThesePos=find(mtmp>(mBL*(sqrt(0.5))));
    % useThesePos=find(mtmp>(min(mtmp)*((10^(thresh_vals(i)/10))*s_fact)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % x=mtmp; % threshold=10; % original- 2.5
    % % Outlier detection
    % [outlierIdx, outlierVals, zScores] = findOutliersMAD(x);
    % useThesePos= outlierIdx;
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
function [outlierIdx, outlierVals, zScores] = findOutliersMAD(x)
% FINDOUTLIERSMAD Detects outliers using Median Absolute Deviation (MAD).
%
%   [outlierIdx, outlierVals, zScores] = findOutliersMAD(x, threshold)
%
%   Inputs:
%       x         - Input data vector (numeric)
%       threshold - Robust z-score cutoff (default = 3.5)
%
%   Outputs:
%       outlierIdx  - Indices of detected outliers
%       outlierVals - Outlier values from x
%       zScores     - Robust z-scores for all data points
%
%   Notes:
%       - Works for skewed/non-normal distributions
%       - Recommended threshold: 3.5 (Iglewicz & Hoaglin, 1993)
%
    % 
    % if nargin < 2
    %     threshold = 2.5; % Default threshold
    % end
    % 
    % Ensure column vector
    x = x(:);
    % Transformation
    % indicator=skewness(x,0);
    % disp(indicator)
    % if indicator>2
    %    x=log(x);
    % end
    % Compute median and raw MAD (no scaling)
    med = median(x);
    mad_raw = mad(x, 1);  % Pure MAD without normality scaling
    
    if mad_raw == 0
        warning('MAD is zero; cannot compute robust z-scores. Data may have no variability.');
        outlierIdx = [];
        outlierVals = [];
        zScores = zeros(size(x));
        return;
    end
    
    % Compute robust z-scores
    zScores = 0.6745 * (x - med) / mad_raw;
    % max_1=max(zScores);
    % max_2=max(zScores(zScores<max(zScores)));
    % max_3=zScores(3);
    % threshold=((max_1+max_2+max_3)/3)*0.8;
    % threshold= 0.5 * max(zScores);
    % threshold=((max(zScores)+min(zScores))/2)*0.7;
    % Identify outliers
    % threshold=2.5;
    % threshold=0.5*power_temp;
    outlierIdx = find((zScores) > threshold);
    outlierVals = x(outlierIdx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mBL = getMeanBaseline(gaborInfo,tOffset,periodS,burstFreqRangeHz)
% SCALE  =1;
FREQ   =2;
POS    =3;
MODULUS=4; %mp2tf uses modulus, although gabor uses amplitude

numTrials = size(gaborInfo,1);
mBLList = [];
for i=1:numTrials
    x = squeeze(gaborInfo(i,:,:));
%   o = x(:,SCALE);
    f = x(:,FREQ);
    t = x(:,POS)+tOffset;
    m = x(:,MODULUS);
    
    goodAtomsInT = intersect(find(t>=periodS(1)),find(t<periodS(2))); % Find Atoms whos time center lies in the stimulus period
    goodAtomsInF = intersect(find(f>=burstFreqRangeHz(1)),find(f<burstFreqRangeHz(2))); % Find Atoms freq center lies in the band
    goodAtomsInTF = intersect(goodAtomsInF,goodAtomsInT);
    
    if ~isempty(goodAtomsInTF)
        mtmp = m(goodAtomsInTF);
        mBLList=cat(2,mBLList,mtmp(mtmp==max(mtmp)));
    end
end
mBL=mean(mBLList);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%