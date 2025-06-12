length_injected=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4];
load('LFP_synth_gamma_all_methods_elec41_ORI_45.mat'); % MP, OMP-GEAR
load('LFP_synth_gamma_rest_methods_elec41_ORI_45.mat'); % Hilbert-Transform, Wavelet-transform, Feingold Algorithm
load('timeVals.mat')
f=figure;
f.WindowState="Maximized";
plotHandles=getPlotHandles(3,1,[0.06 0.08 0.4 0.92],0.04,0.08,0);
% plotHandles_2=getPlotHandles(2,1,[0.52 0.55 0.4 0.43],0.04,0.03,0);
% Raw LFP Trace as well
plotHandles_2=getPlotHandles(3,1,[0.52 0.55 0.42 0.43],0.04,0.05,0);
plotHandles_3=getPlotHandles(1,1,[0.52 0.08 0.42 0.4],0.04,0.08,0);
% Below plot handles for altenrate figure-1 
% plotHandles=getPlotHandles(3,1,[0.06 0.06 0.2 0.92],0.04,0.08,0);
% plotHandles_2=getPlotHandles(3,1,[0.32 0.55 0.6 0.43],0.04,0.03,0);
% plotHandles_3=getPlotHandles(1,1,[0.32 0.06 0.6 0.45],0.04,0.08,0);
% Can include RAW LFP Trace
% Analysis of the performance of the algorithms
% OMP-GEAR- only first column
median_length_OMP_gear_t_1=zeros(1,8);
SEM_length_OMP_gear_t_1=zeros(1,8);
for i=1:8
    % Filter lengths less than 0.8
    length_accumulator_omp_gear{i,1}(length_accumulator_omp_gear{i,1}>0.8)=[];
    median_length_OMP_gear_t_1(i)=median(length_accumulator_omp_gear{i,1});
    SEM_length_OMP_gear_t_1(i)=getSEMedian(length_accumulator_omp_gear{i,1}); 
end
% MP- only first column
median_length_MP_t_1=zeros(1,8);
SEM_length_MP_t_1=zeros(1,8);
for i=1:8
    % Filter lengths less than 0.8
    length_accumulator_MP{i,1}(length_accumulator_MP{i,1}>0.8)=[];
    median_length_MP_t_1(i)=median(length_accumulator_MP{i,1});
    SEM_length_MP_t_1(i)=getSEMedian(length_accumulator_MP{i,1}); 
end
% Hilbert- only first column
median_length_hilbert_t_1=zeros(1,8);
SEM_length_hilbert_t_1=zeros(1,8);
for i=1:8
    % Filter lengths less than 0.8
    length_gatherer_hilbert{i,1}(length_gatherer_hilbert{i,1}>0.8)=[];
    median_length_hilbert_t_1(i)=median(length_gatherer_hilbert{i,1});
    SEM_length_hilbert_t_1(i)=getSEMedian(length_gatherer_hilbert{i,1}); 
end
% Wavelet- only first column
median_length_wavelet_t_1=zeros(1,8);
SEM_length_wavelet_t_1=zeros(1,8);
for i=1:8
    % Filter lengths less than 0.8
    length_gatherer_wavelet{i,1}(length_gatherer_wavelet{i,1}>0.8)=[];
    median_length_wavelet_t_1(i)=median(length_gatherer_wavelet{i,1});
    SEM_length_wavelet_t_1(i)=getSEMedian(length_gatherer_wavelet{i,1}); 
end
% Feingold Algorithm- only first column
% median_length_feingold_t_1=zeros(1,8);
% SEM_length_feingold_t_1=zeros(1,8);
% for i=1:8
%     % Filter lengths less than 0.8
%     length_gatherer_feingold{i,1}(length_gatherer_feingold{i,1}>0.8)=[];
%     median_length_feingold_t_1(i)=median(length_gatherer_feingold{i,1});
%     SEM_length_feingold_t_1(i)=getSEMedian(length_gatherer_feingold{i,1}); 
% end
% Plot the results
subplot(plotHandles(3,1))
errorbar(length_injected,median_length_MP_t_1,SEM_length_MP_t_1,'-o','LineWidth',4,'Color', 'm');
hold on;
% 'm'- MP
errorbar(length_injected,median_length_OMP_gear_t_1,SEM_length_OMP_gear_t_1,'-^','LineWidth',2,'Color', 'g');
% 'g'- OMP-GEAR
hold on;
errorbar(length_injected,median_length_hilbert_t_1,SEM_length_hilbert_t_1,'-s','LineWidth',2,'Color', [0.58, 0.0, 0.83]);
% [0.58, 0.0, 0.83]- Hilbert
hold on;
errorbar(length_injected,median_length_wavelet_t_1,SEM_length_wavelet_t_1,'-d','LineWidth',2,'Color', 'r');
hold on;
% errorbar(length_injected,median_length_feingold_t_1,SEM_length_feingold_t_1,'-x','LineWidth',2,'Color', 'b');
% plot y=x line (plane)- black color
plot([0 0.42],[0 0.42],'--k','LineWidth',0.7);
xlabel('Injected length (s)');
ylabel('Estimated length (s)');
legend('MP','OMP-GEAR','Hilbert','Wavelet');
ylim([0 0.42]);
xlim([0 0.42]);
h3=annotation("textbox",[.004 .22 .1 .1],'String','C','FontSize',24,'FontWeight','Bold','EdgeColor','none','FontName','Helvetica');
h3.FitBoxToText = 'on';
% set gca font size to 18
% set(gca,'FontSize',18);
% Calculate R^2 values
mdl = fitlm(length_injected, median_length_OMP_gear_t_1);
R_squared_omp = mdl.Rsquared.Ordinary;
mdl=fitlm(length_injected, median_length_MP_t_1);
R_squared_mp=mdl.Rsquared.Ordinary;
mdl=fitlm(length_injected, median_length_hilbert_t_1);
R_squared_hilbert=mdl.Rsquared.Ordinary;
mdl=fitlm(length_injected, median_length_wavelet_t_1);
R_squared_wavelet=mdl.Rsquared.Ordinary;
disp(['R^2 value for OMP-GEAR is ',num2str(R_squared_omp)]);
disp(['R^2 value for MP is ',num2str(R_squared_mp)]);
disp(['R^2 value for Hilbert is ',num2str(R_squared_hilbert)]);
disp(['R^2 value for Wavelet is ',num2str(R_squared_wavelet)]);
% generate a power spectral density plot of analogData_accumulator{6,1}, analogData_accumulator{6,2} using chronux toolbox
% for the injected length of 0.3
params.Fs=250;
params.fpass=[0 100];
params.tapers=[1 1];
params.trialave=1;
params.pad=-1;
params.error=[2 0.05];
s1_index_bl=find(timeVals>=0.25,1);
s2_index_bl=find(timeVals>=0.75,1);
select_index=6;
[S,f]=mtspectrumc(analogData_accumulator{select_index,1}(:,s1_index_bl:s2_index_bl)',params);
[S0,f0]=mtspectrumc(analogData_accumulator{select_index,2}(:,s1_index_bl:s2_index_bl)',params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(1,1))
plot(f,(log10(S)),'LineWidth',2,'Color','k');
hold on;
plot(f0,log10(S0),'LineWidth',2,'Color','g');
ylabel('Power (log(\muV^{2}/Hz))');
% set(gca,'FontSize',18);
h2=annotation("textbox",[.004 .54 .1 .1],'String','B','FontSize',24,'FontWeight','Bold','EdgeColor','none','FontName','Helvetica');
h2.FitBoxToText = 'on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(2,1));
s1_index_bl=find(timeVals>=-0.5,1);
s2_index_bl=find(timeVals>=0,1);
[S_bl,f_bl]=mtspectrumc(analogData_accumulator{select_index,1}(:,s1_index_bl:s2_index_bl)',params);
[S0_bl,f0_bl]=mtspectrumc(analogData_accumulator{select_index,2}(:,s1_index_bl:s2_index_bl)',params);
plot(f,10*log10(S./S_bl),'LineWidth',2,'Color','k');
hold on;
plot(f0,10*log10(S0./S0_bl),'LineWidth',2,'Color','g');
ylabel('Power(dB)');
% set(gca,'FontSize',18);
xlabel('Frequency (Hz)');
h1=annotation("textbox",[0.004 0.88 .1 .1],'String','A','FontSize',24,'FontWeight','Bold','EdgeColor','none','FontName','Helvetica');
h1.FitBoxToText = 'on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except plotHandles_3 plotHandles_2 plotHandles f
parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
load(fullfile(parent_file_path,'alpaH_info','parameterCombinations.mat'));
%badtrials file is loaded
load(fullfile(parent_file_path,'alpaH_info','badTrials.mat'));
% selected electrode file is loaded-RF data file
load(fullfile(parent_file_path,'alpaH_info','alpaHMicroelectrodeRFData.mat'));
% LFP_data file is loaded
LFP_data_file=dir(fullfile(parent_file_path,"Decimated_8_LFP_data_alpa_H"));
% LFP_non_decimated_path="C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab\data\alpaH\Microelectrode\210817\GRF_002\segmentedData\LFP\elec41.mat";
% load(LFP_non_decimated_path);
% timeVals_non_decimated_path="C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab\data\alpaH\Microelectrode\210817\GRF_002\segmentedData\LFP\lfpInfo.mat";
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
trial_temp=[parameterCombinations{:,:,:,2,3},parameterCombinations{:,:,:,3,3},parameterCombinations{:,:,:,4,3}];
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
subplot(plotHandles(1,1))
plot(f_stim,log10(S_stim),'LineWidth',2,'Color','r');
hold on;
% draw dotted orange lines at 20 and 32 Hz and dotted blue lines at 36 and 65 Hz
color_orange=[0.9, 0.4, 0.0];
% lighter shade of orange
% color_orange=[1 0.5 0];
color_blue=[0 0 1];
line([slow_gamma_freq(1) slow_gamma_freq(1)], [min(log10(S_stim))-0.5 max(log10(S_stim))+0.5], 'Color', color_blue, 'LineStyle', '--','LineWidth',2);
line([slow_gamma_freq(2) slow_gamma_freq(2)], [min(log10(S_stim))-0.5 max(log10(S_stim))+0.5], 'Color', color_blue, 'LineStyle', '--','LineWidth',2);
line([fast_gamma_freq(1) fast_gamma_freq(1)], [min(log10(S_stim))-0.5 max(log10(S_stim))+0.5], 'Color', color_orange, 'LineStyle', '--','LineWidth',2);
line([fast_gamma_freq(2) fast_gamma_freq(2)], [min(log10(S_stim))-0.5 max(log10(S_stim))+0.5], 'Color', color_orange, 'LineStyle', '--','LineWidth',2);
ylim([min(log10(S_stim))-0.4 max(log10(S_stim))+0.4]);
legend('Injected signal','Spontaneous signal','Original signal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(2,1));
plot(f_stim,10*log10(S_stim./S_bl),'LineWidth',2,'Color','r');
hold on;
line([slow_gamma_freq(1) slow_gamma_freq(1)], [min(10*log10(S_stim./S_bl))-3 max(10*log10(S_stim./S_bl))+3], 'Color', color_blue, 'LineStyle', '--','LineWidth',2);
line([slow_gamma_freq(2) slow_gamma_freq(2)], [min(10*log10(S_stim./S_bl))-3 max(10*log10(S_stim./S_bl))+3], 'Color', color_blue, 'LineStyle', '--','LineWidth',2);
line([fast_gamma_freq(1) fast_gamma_freq(1)], [min(10*log10(S_stim./S_bl))-3 max(10*log10(S_stim./S_bl))+3], 'Color', color_orange, 'LineStyle', '--','LineWidth',2);
line([fast_gamma_freq(2) fast_gamma_freq(2)], [min(10*log10(S_stim./S_bl))-3 max(10*log10(S_stim./S_bl))+3], 'Color', color_orange, 'LineStyle', '--','LineWidth',2);
ylim([min(10*log10(S_stim./S_bl))-2.8 max(10*log10(S_stim./S_bl))+2.8]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend('Injected signal','Spontaneous signal','Original signal');
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
    % LFP_signal_non_deci=analogData(trial_temp,:);
    % Plot the original LFP data signal
    LFP_signal_one_trial=LFP_signal((trial),:);
    % LFP_signal_one_trial_non_deci=LFP_signal((trial),:);
    % filter between 20 Hz and 32 Hz
    slow_gamma_freq_burst_range=[24 29];
    [b,a]=butter(4,slow_gamma_freq_burst_range/(Fs/2),'bandpass');
    LFP_signal_one_trial_sg=filtfilt(b,a,LFP_signal_one_trial);
    % filter between 36 Hz and 65 Hz
    fast_gamma_freq_burst_range=[38 54]; % Optimize the frequency range to localize the bursts
    [b,a]=butter(4,fast_gamma_freq_burst_range/(Fs/2),'bandpass');
    LFP_signal_one_trial_fg=filtfilt(b,a,LFP_signal_one_trial);
    % plot the signal
    %%%%%%%%%%%%%%%%%%%%%%% SLOW GAMMA %%%%%%%%%%%%%%%%%%%
    % load(timeVals_non_decimated_path);
    subplot(plotHandles_2(3,1));
    % subplot(plotHandles_2(2,1));
    plot(timeVals,LFP_signal_one_trial_sg,'-b');
    ylim([min(LFP_signal_one_trial_sg)+5,max(LFP_signal_one_trial_sg)+5])
    % Keepeing the y-axis same (of both fast gamma and slow gamma) for better comparison
    % Fast gamma has a higher amplitude than slow gamma (for this particular trial)
    % ylim([min(LFP_signal_one_trial_fg)-5,max(LFP_signal_one_trial_fg)+5])
    hold on;
    % set(gca,'FontSize',18);
    %plot(timeVals,LFP_signal_one_trial,'-k');

    % title('Slow gamma- sample trial- Bursts estimation by MP');
    %%%%%%%%%%%%%%%%%%%%% FAST GAMMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles_2(2,1));
    % subplot(plotHandles_2(1,1));
    plot(timeVals,LFP_signal_one_trial_fg,'-','Color',color_orange);
    % ylim([min(LFP_signal_one_trial_fg)-5,max(LFP_signal_one_trial_fg)+5])
    ylim([min(LFP_signal_one_trial_sg)+5,max(LFP_signal_one_trial_sg)+5])
    hold on;
    ylabel("Amplitude(\muV)")
    % set(gca,'FontSize',18);
    % title('Fast gamma- sample trial- Bursts estimation by MP');
    %%%%%%%%%%%%%%%%%%%% LFP signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles_2(1,1));
    % plot(timeVals,LFP_signal_one_trial,'-k');
    % load(timeVals_non_decimated_path);
    plot(timeVals,LFP_signal_one_trial,'-k','LineWidth',1.5);
    ylim([min(LFP_signal_one_trial)-10,max(LFP_signal_one_trial)+10])
    xlim([0 1]);
    hold on;
    line([0.25 0.25], [-420 220], 'Color', 'm', 'LineStyle', '--');
    line([0.75 0.75], [-420 220], 'Color', 'm', 'LineStyle', '--');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    thresholdFraction=0.5;
    num_iterations=120;
    dict_Size=2500000;
    jj=selected_elec_pos;
    displayFlag=1;
    % current elec- 41 (pos-37)
    data_temp=LFP_signal_one_trial;
    gabor_SF=gaborInfo_accumulator_alpaH{current_SF,current_ORI,jj};
    gabor_SF_one_trial=gabor_SF(trial,:,:);
    header_SF=header_accumulator_alpaH{current_SF,current_ORI,jj};
    header_SF_one_trial=header_SF(trial,:,:);
    %Slow gamma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma_freq=[20 65];
    diffPower=getChangeInPower_single(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,gamma_freq);
    thresholdFactor=sqrt(thresholdFraction*diffPower);
    [length_measured_broad,freq_measured_broad,time_measured_broad,~,~,~]= getBurstLengthMP_2(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,gamma_freq,num_iterations,0.9,dict_Size,gabor_SF_one_trial,header_SF_one_trial,plotHandles_3);
    subplot(plotHandles_3(1,1));
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
    displayFlag=0;
    % Slow gamma burst computation
    data_temp=LFP_signal;
    diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,slow_gamma_freq);
    thresholdFactor=sqrt(thresholdFraction*diffPower);
    [length_measured_sg,freq_measured_sg,time_measured_sg,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,num_iterations,0.9,dict_Size,gabor_SF,header_SF);
    % Fast gamma burst computation
    diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
    thresholdFactor=sqrt(thresholdFraction*diffPower);
    [length_measured_fg,freq_measured_fg,time_measured_fg,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,num_iterations,0.9,dict_Size,gabor_SF,header_SF);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles_2(2,1));
    hold on;
    plot_burst_lengths(length_measured_fg,time_measured_fg,trial,'m');
    subplot(plotHandles_2(3,1));
    hold on;
    plot_burst_lengths(length_measured_sg,time_measured_sg,trial,'g');
    h4=annotation("textbox",[0.47 0.88 .1 .1],'String','D','FontSize',24,'FontWeight','Bold','EdgeColor','none','FontName','Helvetica');
    h4.FitBoxToText = 'on';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles_3(1,1));
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
    % set(gca,'FontSize',18);
    c=colorbar;
    c.Position = c.Position + [0.06 0 0.005 0];
    c.Ticks  = [-3 -2 -1 0 1 2 3];
    c.FontSize = 10;
    % c.Label.FontSize = 14;
    box(c,'off');
    h5=annotation("textbox",[0.47 0.42 .1 .1],'String','E','FontSize',24,'FontWeight','Bold','EdgeColor','none','FontName','Helvetica');
    h5.FitBoxToText = 'on';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %title(['MP- Sample Trial - ',num2str(trial)])
    axis([0 1 0 85]);
    xlabel('Time(s)')
    ylabel('Frequency(Hz)')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set_axis_ticks_fontsize(plotHandles,22,18,1);
    set_axis_ticks_fontsize(plotHandles,22,18,2);
    set_axis_ticks_fontsize(plotHandles,22,18,3);
    set_axis_ticks_fontsize(plotHandles_2,22,18,1);
    set_axis_ticks_fontsize(plotHandles_2,22,18,2);
    set_axis_ticks_fontsize(plotHandles_2,22,18,3);
    set_axis_ticks_fontsize(plotHandles_3,22,18,1);
    t1=annotation('textbox',...
    [0.97 0.55 0.06 0.06],...
    'String',{'Slow \gamma'},...
    'Rotation',90,...
    'FontWeight','bold',...
    'FontSize',16,...
    'FontName','Helvetica',...
    'EdgeColor',[1 1 1]);
    t1.FitBoxToText = 'on';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t2=annotation('textbox',...
    [0.97 0.71 0.06 0.06],...
    'String',{'Fast \gamma'},...
    'Rotation',90,...
    'FontWeight','bold',...
    'FontSize',16,...
    'FontName','Helvetica',...
    'EdgeColor',[1 1 1]);
     t2.FitBoxToText = 'on';