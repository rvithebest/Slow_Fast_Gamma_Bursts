clc;clear all;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
load(fullfile(parent_file_path,'alpaH_info','parameterCombinations.mat'));
%badtrials file is loaded
load(fullfile(parent_file_path,'alpaH_info','badTrials.mat'));
% selected electrode file is loaded-RF data file
load(fullfile(parent_file_path,'alpaH_info','alpaHMicroelectrodeRFData.mat'));
% LFP_data file is loaded
LFP_data_file=dir(fullfile(parent_file_path,"Decimated_8_LFP_data_alpa_H"));
load('gamma_duration_alpaH_MP.mat');
% natrisort the LFP data-remain in form of struct
LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'})); %remove . and ..
%I want all the files in the folder to be sorted by their name in order 1,2,3... (not 1,10,11,..)
LFP_data_file = natsortfiles({LFP_data_file.name});
stimulusPeriodS=[0.25 0.75];
baselinePeriodS=[-0.5 0];
%%%%%%%% Used for calculating the change in power %%%%%%%%%%%%%%
fast_gamma_freq=[36 65];
%%%%%%% for alpaH - slow gamma range is taken from 25 to 40 Hz%%%%%%%%%
slow_gamma_freq=[20 32];
selected_elec=highRMSElectrodes;
selected_elec_LFP=selected_elec(1:77);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%j=1;% counter for electrode position
electrode_num=length(selected_elec_LFP);
%%%%%%%%%%%%%%%%%%%%%%%%%%% LFP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 32- possible elec
selected_elec_pos=37;
i=selected_elec_LFP(selected_elec_pos);
load(fullfile('Decimated_8_LFP_Data_alpa_H',LFP_data_file{i}));%analogDataDecimatedDecimated
load('timeVals.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
current_SF=2; % 1 CPD
current_ORI=3; % 45 deg
trial=1; % first trial of SF- 1 cpd and 45 deg
% SF- 1 cpd and 157.5 deg- trial num=6, 13, 17
%%%%%%%%%%%%%% BAND-PASS FILTERING %%%%%%%%%%%%%%%%%%%%%%
trial_temp=[parameterCombinations{:,:,:,current_SF,current_ORI}];
trial_temp=setdiff(trial_temp,badTrials);
LFP_signal=analogDataDecimated(trial_temp,:);
% Plot the original LFP data signal
LFP_signal_one_trial=LFP_signal((trial),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresholdFraction=0.5;
dict_Size=2500000;
jj=selected_elec_pos;
displayFlag=0;
% current elec- 41 (pos-37)
data_temp=LFP_signal_one_trial;
%%%%%%%%%%%%%%%%%% Slow gamma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_freq=[20 65];
num_iterations=600;
diffPower=getChangeInPower_single(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,gamma_freq);
thresholdFactor=sqrt(thresholdFraction*diffPower);
[~,~,~,gabor_temp_one_trial,header_temp_one_trial,~]= getBurstLengthMP_non_filter(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,gamma_freq,num_iterations,0.9,dict_Size,[],[]);
save('MP_non_filt_single_trial_results.mat','gabor_temp_one_trial','header_temp_one_trial');