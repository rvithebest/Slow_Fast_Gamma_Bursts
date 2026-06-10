f=figure;
f.WindowState="Maximized";
plotHandles_3=getPlotHandles(3,1,[0.48 0.55 0.46 0.39],0.04,0.05,0);
plotHandles_4=getPlotHandles(1,1,[0.48 0.08 0.46 0.36],0.04,0.08,0);
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
trial_temp=[parameterCombinations{:,:,:,2,9}];
trial_temp=setdiff(trial_temp,badTrials);
LFP_signal=analogDataDecimated(trial_temp,:);
Fs=250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PSD plot of LFP signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.Fs=Fs;
params.fpass=[0 100];
params.tapers=[1 1];
params.trialave=1;
params.pad=-1;
params.err=[2 0.05];
s1_index=find(timeVals>=0.25,1);
s2_index=find(timeVals>=0.75,1);
LFP_signal_stim=LFP_signal(:,s1_index:s2_index);
[S_stim,f_stim]=mtspectrumc(LFP_signal_stim',params);
s1_index_bl=find(timeVals>=-0.5,1);
s2_index_bl=find(timeVals>=0,1);
LFP_signal_bl=LFP_signal(:,s1_index_bl:s2_index_bl);
[S_bl,f_bl]=mtspectrumc(LFP_signal_bl',params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_orange=[0.9, 0.4, 0.0];
color_blue=[0 0 1];
current_SF=2; % 1 CPD
current_ORI=3; % 45 deg
trial=1; % first trial of SF- 1 cpd and 45 deg
% SF- 1 cpd and 157.5 deg- trial num=6, 13, 17
%%%%%%%%%%%%%% BAND-PASS FILTERING %%%%%%%%%%%%%%%%%%%%%%
trial_temp=[parameterCombinations{:,:,:,current_SF,current_ORI}];
trial_temp=setdiff(trial_temp,badTrials);
LFP_signal=analogDataDecimated(trial_temp,:);
% LFP data signal
LFP_signal_one_trial=LFP_signal((trial),:);
% Filtered version of the LFP signal in the slow gamma and fast gamma range
slow_gamma_freq_burst_range=[24 29];
[b,a]=butter(4,slow_gamma_freq_burst_range/(Fs/2),'bandpass');
LFP_signal_one_trial_sg=filtfilt(b,a,LFP_signal_one_trial);
fast_gamma_freq_burst_range=[38 54]; % Optimize the frequency range to localize the bursts
[b,a]=butter(4,fast_gamma_freq_burst_range/(Fs/2),'bandpass');
LFP_signal_one_trial_fg=filtfilt(b,a,LFP_signal_one_trial);
%%%%%%%%%%%%%%%%%%%%%%% SLOW GAMMA %%%%%%%%%%%%%%%%%%%
subplot(plotHandles_3(3,1));
plot(timeVals,LFP_signal_one_trial_sg,'-b');
ylim([min(LFP_signal_one_trial_sg)+5,max(LFP_signal_one_trial_sg)+5]);
%%%%%%%%%%%%%%%%%%%%% FAST GAMMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_3(2,1));
plot(timeVals,LFP_signal_one_trial_fg,'-','Color',color_orange);
ylim([min(LFP_signal_one_trial_sg)+5,max(LFP_signal_one_trial_sg)+5]);
hold on;
ylabel("Amplitude(\muV)")
%%%%%%%%%%%%%%%%%%%% LFP signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_3(1,1));
plot(timeVals,LFP_signal_one_trial,'-k','LineWidth',1.5);
ylim([min(LFP_signal_one_trial)-10,max(LFP_signal_one_trial)+10])
xlim([0 1]);
hold on;
line([0.25 0.25], [-520 220], 'Color', 'm', 'LineStyle', '--','LineWidth',2);
line([0.75 0.75], [-520 220], 'Color', 'm', 'LineStyle', '--','LineWidth',2);
line([-0.5 -0.5], [-520 220], 'Color',[0.23,0.51,0], 'LineStyle', '--','LineWidth',2);
line([0 0], [-520 220], 'Color',[0.23,0.51,0], 'LineStyle', '--','LineWidth',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresholdFraction=0.5;
num_iterations=120;
dict_Size=2500000;
jj=selected_elec_pos;
displayFlag=1;
data_temp=LFP_signal_one_trial;
gabor_SF=gaborInfo_accumulator_alpaH{current_SF,current_ORI,jj};
gabor_SF_one_trial=gabor_SF(trial,:,:);
header_SF=header_accumulator_alpaH{current_SF,current_ORI,jj};
header_SF_one_trial=header_SF(trial,:,:);
%%%%%%%%%%%%%%%%%% Slow gamma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_freq=[20 65];
diffPower=getChangeInPower_single(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,gamma_freq);
thresholdFactor=sqrt(thresholdFraction*diffPower);
idx_1=1;idx_2=1;
[~,~,~,~,~,~]= getBurstLengthMP_plot(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,gamma_freq,num_iterations,0.9,dict_Size,gabor_SF_one_trial,header_SF_one_trial,plotHandles_4,idx_1,idx_2);
subplot(plotHandles_4(1,1));
hold on;
% set y-axis to [0 100] and x axis to the timeVals
axis([timeVals(1) timeVals(end) 0 100]);
% plot the  horizontal line (Across time Vals)- frequency ranges using dashed white lines(sloww gamma) and dashed black lines(fast gamma)
line([timeVals(1) timeVals(end)],[slow_gamma_freq(1) slow_gamma_freq(1)],'Color','k','LineStyle','--','LineWidth',2);
line([timeVals(1) timeVals(end)],[slow_gamma_freq(2) slow_gamma_freq(2)],'Color','k','LineStyle','--','LineWidth',2);
% plot the  horizontal line (Across time Vals)- frequency ranges using dashed white lines(sloww gamma) and dashed black lines(fast gamma)
line([timeVals(1) timeVals(end)],[fast_gamma_freq(1) fast_gamma_freq(1)],'Color','k','LineStyle','--','LineWidth',2);
line([timeVals(1) timeVals(end)],[fast_gamma_freq(2) fast_gamma_freq(2)],'Color','k','LineStyle','--','LineWidth',2);
% plot the burst length at time burst center and y-axis frequency_measured
load("Detected_burst_data_MP_sample_elec.mat","length_measured_sg","freq_measured_sg","time_measured_sg","length_measured_fg","freq_measured_fg","time_measured_fg");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_3(2,1));
hold on;
plot_burst_lengths(length_measured_fg,time_measured_fg,trial,'m');
subplot(plotHandles_3(3,1));
hold on;
plot_burst_lengths(length_measured_sg,time_measured_sg,trial,'g');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_4(1,1));
hold on;
for ii = 1:length(length_measured_sg{trial})
    if length_measured_sg{trial}(ii)>0.8
        continue
    end
    start_time = time_measured_sg{trial}(ii) - (length_measured_sg{trial}(ii)/2);
    end_time = time_measured_sg{trial}(ii) + (length_measured_sg{trial}(ii)/2);
    % plot the burst length at time burst center- same as before and y-axis frequency_measured
    plot([start_time,end_time],[freq_measured_sg{trial}(ii),freq_measured_sg{trial}(ii)],'Color','k','LineWidth',2);
    % plot the burst length at time burst center- same as before and y-axis frequency_measured
    plot(time_measured_sg{trial}(ii),freq_measured_sg{trial}(ii),'ko','MarkerSize',5,'MarkerFaceColor','g');
end
for ii=1:length(length_measured_fg{trial})
    if length_measured_fg{trial}(ii)>0.8
        continue
    end
    start_time = time_measured_fg{trial}(ii) - length_measured_fg{trial}(ii)/2;
    end_time = time_measured_fg{trial}(ii) + length_measured_fg{trial}(ii)/2;
    % plot the burst length at time burst center- same as before and y-axis frequency_measured
    plot([start_time,end_time],[freq_measured_fg{trial}(ii),freq_measured_fg{trial}(ii)],'Color','k','LineWidth',2);
    % plot the burst length at time burst center- same as before and y-axis frequency_measured
    plot(time_measured_fg{trial}(ii),freq_measured_fg{trial}(ii),'ko','MarkerSize',5,'MarkerFaceColor','m');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % draw a white dotted line at 0.25 and 0.75
line([0.25 0.25], [0 100], 'Color','m', 'LineStyle', '--','LineWidth',2);
line([0.75 0.75], [0 100], 'Color', 'm', 'LineStyle', '--','LineWidth',2);
line([-0.5 -0.5], [-100 100], 'Color',[0.23,0.51,0], 'LineStyle', '--','LineWidth',2);
line([0 0], [-100 100], 'Color',[0.23,0.51,0], 'LineStyle', '--','LineWidth',2);
c=colorbar;
c.Position = c.Position + [0.06 0 0.005 0];
c.Ticks  = [-4 -2 0 2 4];
c.FontSize = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis([-0.1 1 0 80]);
clim([-4 4]);
xlabel('Time(s)');
ylabel('Frequency(Hz)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set_axis_ticks_fontsize(plotHandles_3,22,18,1);
set_axis_ticks_fontsize(plotHandles_3,22,18,2);
set_axis_ticks_fontsize(plotHandles_3,22,18,3);
set_axis_ticks_fontsize(plotHandles_4,22,18,1);
t1=annotation('textbox',...
[0.97 0.53 0.06 0.06],...
'String',{'Slow \gamma'},...
'Rotation',90,...
'FontWeight','bold',...
'FontSize',16,...
'Color','b',...
'FontName','Helvetica',...
'EdgeColor',[1 1 1]);
t1.FitBoxToText = 'on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t2=annotation('textbox',...
[0.97 0.69 0.06 0.06],...
'String',{'Fast \gamma'},...
'Rotation',90,...
'FontWeight','bold',...
'FontSize',16,...
'Color',color_orange,...
'FontName','Helvetica',...
'EdgeColor',[1 1 1]);
 t2.FitBoxToText = 'on';
labels = {'A','B'};
x_positions = [0.435];
y_positions = [0.87,0.38];  
k = 1;
for j = 1:length(y_positions)
    for i = 1:length(x_positions)
        annotation('textbox', ...
            [x_positions(i), y_positions(j), 0.1, 0.1], ...
            'String', labels{k}, ...
            'FontSize', 28, ...
            'FontWeight', 'Bold', ...
            'EdgeColor', 'none', ...
            'FontName', 'Helvetica');
        k = k + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R_squared_val, RMSE_val] = compute_R_squared_val(x,y)
    % Computes the R^2 value
    x=x'; y=y';
    SSresid = sum((y - x).^2);        % residual sum of squares
    SStotal = sum((y - mean(y)).^2);  % total variance of y
    R_squared_val = 1 - SSresid / SStotal;
    % RMSE value
    RMSE_val = sqrt(mean((y - x).^2));
end