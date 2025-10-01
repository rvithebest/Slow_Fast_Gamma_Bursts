load('Model_FG_3/simulationResults.mat');
t = JSpop.EIpairs.t;
rE = JSpop.EIpairs.R(1:end/2,:);
rI = JSpop.EIpairs.R(end/2+1:end,:);
t=t(1:end-1);
lfp=lfp(:,1:end-1);
t_sel_analysis = t>=0.25 & t<=1;
lfp_fft_trialwise = abs(fft(lfp(:,t_sel_analysis),[],2));
Fs=1/(t(2)-t(1));
freq_axis=0:Fs/length(t(t_sel_analysis)):Fs*(1-1/length(t(t_sel_analysis)));
mean_lfp_fft = mean(lfp_fft_trialwise, 1);
decimationFactor=80;
lfp_decimated = zeros(size(lfp,1), size(lfp,2)/decimationFactor);
for numTrial=1:size(lfp,1)
        lfp_decimated(numTrial,:)=resample(lfp(numTrial,:)',1,decimationFactor)'; 
end
timeVals_decimated = downsample(t', decimationFactor)';
%%%%%%%%%%%%%%% MP length analysis %%%%%%%%%%%%%%%%%%
displayFlag=0;
stimulusPeriodS=[0.25 1];
baselinePeriodS=[-0.75 0];
gamma_freq=[40 70];
thresholdFraction=0.5;
num_iterations=120;
dict_size=2500000;
diffPower=getChangeInPower(lfp_decimated,timeVals_decimated,stimulusPeriodS,baselinePeriodS,gamma_freq);
thresholdFactor=sqrt(thresholdFraction*diffPower);
[length_temp_sg,~,~,gabor_temp,header_temp,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
length_temp_all_trials=[];
for ii=1:length(length_temp_sg)
    length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}'];
end