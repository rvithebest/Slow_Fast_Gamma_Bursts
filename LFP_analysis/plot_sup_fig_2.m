parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
load(fullfile(parent_file_path,'alpaH_info','parameterCombinations.mat'));
%badtrials file is loaded
load(fullfile(parent_file_path,'alpaH_info','badTrials.mat'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('LFP_synth_gamma_all_methods_elec41_all_ORI.mat')
load('timeVals.mat');
f=figure;
f.WindowState="Maximized";
plotHandles_1=getPlotHandles(1,1,[0.06 0.08 0.35 0.55],0.04,0.08,0);
plotHandles_2=getPlotHandles(2,1,[0.06 0.7 0.35 0.28],0.04,0.01,0);
plotHandles_3=getPlotHandles(3,1,[0.48 0.08 0.46 0.9],0.04,0.07,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slow_gamma_freq=[20 36];
fast_gamma_freq=[36 65];
displayFlag=0;
stimulusPeriodS=[0.25 0.75];
baselinePeriodS=[-0.5 0];
thresholdFraction=0.5;
cutoff_criteria=0.04;
trial_idx_gatherer=cell(1,9);
% ORI-wise burst measurement
for ORI=1:8
    trial_all_ORI=(parameterCombinations{:,:,:,2,9}); % SF-1 cpd, all ORI
    trial_all_ORI=setdiff(trial_all_ORI,badTrials);
    trial_ORI_temp=(parameterCombinations{:,:,:,2,ORI});
    trial_ORI_temp=setdiff(trial_ORI_temp,badTrials);
    % Find those indices and store in trial_idx_temp
    trial_idx_temp=find(ismember(trial_all_ORI,trial_ORI_temp));
    trial_idx_gatherer{ORI}=trial_idx_temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
median_length_OMP_gear_t_1=zeros(1,8);
SEM_length_OMP_gear_t_1=zeros(1,8);
recall_OMP_gear_t_1=zeros(1,8);
precision_OMP_gear_t_1=zeros(1,8);
F1_score_OMP_gear_t_1=zeros(1,8);
for i=1:8
    length_temp_all=[];
    freq_temp_all=[];
    for ORI=1:8
        data_temp=analogData_accumulator{i}(trial_idx_gatherer{ORI},:);
        freq_burst_temp=freqBurst_accumulator{i}(trial_idx_gatherer{ORI});
        gabor_temp=gabor_accmulator_omp_gear{i}(trial_idx_gatherer{ORI},:,:);
        % Get burst durations
        % Slow gamma
        diffPower=getChangeInPower(data_temp,timeVals_accumulator{i},stimulusPeriodS,baselinePeriodS,slow_gamma_freq);
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        [length_sg,freq_sg,time_sg,~,~,~]=getBurstLength_all(data_temp,timeVals_accumulator{i},thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,120,0.9,1500000,gabor_temp,[],'OMP-GEAR');
        for ii=1:length(length_sg)
            if isempty(length_sg{ii})
                continue;
            end
            % reject idx longer than 0.8 seconds
            reject_idx=length_sg{ii}>0.8;
            length_sg{ii}(reject_idx)=[];
            freq_sg{ii}(reject_idx)=[];
            length_temp_all=[length_temp_all, length_sg{ii}'];
            freq_temp_all=[freq_temp_all, freq_sg{ii}'];
        end
        % Fast gamma
        diffPower=getChangeInPower(data_temp,timeVals_accumulator{i},stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        [length_fg,freq_fg,time_fg,~,~,~]=getBurstLength_all(data_temp,timeVals_accumulator{i},thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,120,0.9,1500000,gabor_temp,[],'OMP-GEAR');
        for ii=1:length(length_fg)
            if isempty(length_fg{ii})
                continue
            end
            % reject idx longer than 0.8 seconds
            reject_idx=length_fg{ii}>0.8;
            length_fg{ii}(reject_idx)=[];
            freq_fg{ii}(reject_idx)=[];
            length_temp_all=[length_temp_all, length_fg{ii}'];
            freq_temp_all=[freq_temp_all, freq_fg{ii}'];
        end
    end
    % Filter lengths less than 0.8 seconds
    median_length_OMP_gear_t_1(i)=median(length_temp_all);
    SEM_length_OMP_gear_t_1(i)=getSEMedian(length_temp_all,1000);
    total_burst_inj=sum(cell2mat(numBurst_accumulator{i,1}));
    actual_length=length_injected(i);
    corr_bursts=length(length_temp_all(length_temp_all>=(actual_length-cutoff_criteria) & length_temp_all<=(actual_length+cutoff_criteria)));
    recall_OMP_gear_t_1(i)=corr_bursts/total_burst_inj;
    false_bursts=length(length_temp_all)-corr_bursts;
    precision_OMP_gear_t_1(i)=corr_bursts/(corr_bursts+false_bursts);
    F1_score_OMP_gear_t_1(i)=2*(precision_OMP_gear_t_1(i)*recall_OMP_gear_t_1(i))/(precision_OMP_gear_t_1(i)+recall_OMP_gear_t_1(i));
end
% MP
median_length_MP_t_1=zeros(1,8);
SEM_length_MP_t_1=zeros(1,8);
recall_MP_t_1=zeros(1,8);
precision_MP_t_1=zeros(1,8);
F1_score_MP_t_1=zeros(1,8);
for i=1:8
    length_temp_all=[];
    freq_temp_all=[];
    for ORI=1:8
        data_temp=analogData_accumulator{i}(trial_idx_gatherer{ORI},:);
        freq_burst_temp=freqBurst_accumulator{i}(trial_idx_gatherer{ORI});
        gabor_temp=gabor_accmulator_MP{i}(trial_idx_gatherer{ORI},:,:);
        % Get burst durations
        % Slow gamma
        diffPower=getChangeInPower(data_temp,timeVals_accumulator{i},stimulusPeriodS,baselinePeriodS,slow_gamma_freq);
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        [length_sg,freq_sg,time_sg,~,~,~]=getBurstLength_all(data_temp,timeVals_accumulator{i},thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,120,0.9,2500000,gabor_temp,[],'MP');
        for ii=1:length(length_sg)
            if isempty(length_sg{ii})
                continue
            end
            % reject idx longer than 0.8 seconds
            reject_idx=length_sg{ii}>0.8;
            length_sg{ii}(reject_idx)=[];
            freq_sg{ii}(reject_idx)=[];
            length_temp_all=[length_temp_all, length_sg{ii}'];
            freq_temp_all=[freq_temp_all, freq_sg{ii}'];
        end
        % Fast gamma
        diffPower=getChangeInPower(data_temp,timeVals_accumulator{i},stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        [length_fg,freq_fg,time_fg,~,~,~]=getBurstLength_all(data_temp,timeVals_accumulator{i},thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,120,0.9,2500000,gabor_temp,[],'MP');
        for ii=1:length(length_fg)
            if isempty(length_fg{ii})
                continue
            end
            % reject idx longer than 0.8 seconds
            reject_idx=length_fg{ii}>0.8;
            length_fg{ii}(reject_idx)=[];
            freq_fg{ii}(reject_idx)=[];
            length_temp_all=[length_temp_all, length_fg{ii}'];
            freq_temp_all=[freq_temp_all, freq_fg{ii}'];
        end
    end
    if i==6 % For 3oo ms
        subplot(plotHandles_3(2,1));
        hold on;
        bin_size=4; % from 20 to 65 Hz
        % find corresponding bin for each frequency and thus corresponding length
        freq_bins=20:bin_size:65;
        freq_bin_centers=freq_bins(1:end-1)+bin_size/2;
        freq_bins_num=discretize(freq_temp_all,freq_bins);
        % find median length and SEM corresponding to each bin
        median_length_freq_bins=zeros(1,length(freq_bin_centers));
        sem_length_freq_bins=zeros(1,length(freq_bin_centers));
        for k=1:(length(freq_bins)-1)
            temp_indices=find(freq_bins_num==k);
            temp_lengths=length_temp_all(temp_indices);
            median_length_freq_bins(k)=median(temp_lengths);
            sem_length_freq_bins(k)=getSEMedian(temp_lengths,1000);
        end
        errorbar(freq_bin_centers,median_length_freq_bins,sem_length_freq_bins,'-o','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r','Color','r','LineWidth',1.5);
        % plot a horizontal dashed line at y=0.3 seconds
        plot([freq_bins(1) freq_bins(end)], [0.3 0.3], '--k', 'LineWidth', 1);
        xlabel('Frequency (Hz)');
        ylabel('Burst length (s)');
      end
        median_length_MP_t_1(i)=median(length_temp_all);
        SEM_length_MP_t_1(i)=getSEMedian(length_temp_all,1000);
        total_burst_inj=sum(cell2mat(numBurst_accumulator{i,1}));
        actual_length=length_injected(i);
        corr_bursts=length(length_temp_all(length_temp_all>=(actual_length-cutoff_criteria) & length_temp_all<=(actual_length+cutoff_criteria)));
        recall_MP_t_1(i)=corr_bursts/total_burst_inj;
        false_bursts=length(length_temp_all)-corr_bursts;
        precision_MP_t_1(i)=corr_bursts/(corr_bursts+false_bursts);
        F1_score_MP_t_1(i)=2*(precision_MP_t_1(i)*recall_MP_t_1(i))/(precision_MP_t_1(i)+recall_MP_t_1(i));
end
% Hilbert
median_length_hilbert_t_1=zeros(1,8);
SEM_length_hilbert_t_1=zeros(1,8);
recall_hilbert_t_1=zeros(1,8);
precision_hilbert_t_1=zeros(1,8);
F1_score_hilbert_t_1=zeros(1,8);
filter_order=4;
for i=1:8
   length_temp_all=[];
    for ORI=1:8
        data_temp=analogData_accumulator{i}(trial_idx_gatherer{ORI},:);
        freq_burst_temp=freqBurst_accumulator{i}(trial_idx_gatherer{ORI});
        % Get burst durations
        % Slow gamma
        diffPower=getChangeInPower(data_temp,timeVals_accumulator{i},stimulusPeriodS,baselinePeriodS,slow_gamma_freq);
        thresholdFactor=(diffPower*thresholdFraction);
        [length_sg]=getBurstLengthHilbert(data_temp,timeVals_accumulator{i},thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,filter_order);
        for ii=1:length(length_sg)
            if isempty(length_sg{ii})
                continue
            end
            % reject idx longer than 0.8 seconds
            reject_idx=length_sg{ii}>0.8;
            length_sg{ii}(reject_idx)=[];
            length_temp_all=[length_temp_all, length_sg{ii}];
        end
        % Fast gamma
        diffPower=getChangeInPower(data_temp,timeVals_accumulator{i},stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
        thresholdFactor=(diffPower*thresholdFraction);
        [length_fg]=getBurstLengthHilbert(data_temp,timeVals_accumulator{i},thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,filter_order);
        for ii=1:length(length_fg)
            if isempty(length_fg{ii})
                continue
            end
            % reject idx longer than 0.8 seconds
            reject_idx=length_fg{ii}>0.8;
            length_fg{ii}(reject_idx)=[];
            length_temp_all=[length_temp_all, length_fg{ii}];
        end
    end
    median_length_hilbert_t_1(i)=median(length_temp_all);
    SEM_length_hilbert_t_1(i)=getSEMedian(length_temp_all,1000);
    total_burst_inj=sum(cell2mat(numBurst_accumulator{i,1}));
    actual_length=length_injected(i);
    corr_bursts=length(length_temp_all(length_temp_all>=(actual_length-cutoff_criteria) & length_temp_all<=(actual_length+cutoff_criteria)));
    recall_hilbert_t_1(i)=corr_bursts/total_burst_inj;
    false_bursts=length(length_temp_all)-corr_bursts;
    precision_hilbert_t_1(i)=corr_bursts/(corr_bursts+false_bursts);
    F1_score_hilbert_t_1(i)=2*(precision_hilbert_t_1(i)*recall_hilbert_t_1(i))/(precision_hilbert_t_1(i)+recall_hilbert_t_1(i));
end
% Plot the results
subplot(plotHandles_1(1,1))
hold on;
% 'm'- MP
errorbar(length_injected,median_length_OMP_gear_t_1,SEM_length_OMP_gear_t_1,'-^','LineWidth',2,'Color', 'b');
% 'g'- OMP-GEAR
hold on;
errorbar(length_injected,median_length_hilbert_t_1,SEM_length_hilbert_t_1,'-s','LineWidth',2,'Color', [0.58, 0.0, 0.83]);
% [0.58, 0.0, 0.83]- Hilbert
hold on;
errorbar(length_injected,median_length_MP_t_1,SEM_length_MP_t_1,'-o','LineWidth',2,'Color', 'r');
hold on;
% plot y=x line (plane)- black color
plot([0 0.42],[0 0.42],'--k','LineWidth',0.7);
xlabel('Injected length (s)');
ylabel('Estimated length (s)');
legend('OMP-GEAR','Hilbert','MP');
ylim([0 0.42]);
xlim([0 0.42]);
% Calculate R^2 values
[R_squared_OMP, RMSE_OMP] = compute_R_squared_val(length_injected, median_length_OMP_gear_t_1);
[R_squared_MP, RMSE_MP] = compute_R_squared_val(length_injected, median_length_MP_t_1);
[R_squared_hilbert, RMSE_hilbert] = compute_R_squared_val(length_injected, median_length_hilbert_t_1);
disp(['R^2 value for OMP-GEAR is ',num2str(R_squared_OMP)]);
disp(['R^2 value for MP is ',num2str(R_squared_MP)]);
disp(['R^2 value for Hilbert is ',num2str(R_squared_hilbert)]);
disp(['RMSE value for OMP-GEAR is ',num2str(RMSE_OMP)]);
disp(['RMSE value for MP is ',num2str(RMSE_MP)]);
disp(['RMSE value for Hilbert is ',num2str(RMSE_hilbert)]);
subplot(plotHandles_3(3,1));
plot(length_injected, F1_score_OMP_gear_t_1, '-^', 'LineWidth', 2, 'Color', 'b');
hold on;
plot(length_injected, F1_score_hilbert_t_1, '-s', 'LineWidth', 2, 'Color', [0.58, 0.0, 0.83]);
hold on;
plot(length_injected, F1_score_MP_t_1, '-o', 'LineWidth', 2, 'Color', 'r');
xlabel('Injected length (s)');
ylabel('F1-Score');
legend('OMP-GEAR','Hilbert','MP');
% for the injected length of 300 ms
params.Fs=250;
params.fpass=[0 100];
params.tapers=[1 1];
params.trialave=1;
params.pad=-1;
params.error=[2 0.05];
s1_index_stim=find(timeVals>=0.25,1);
s2_index_stim=find(timeVals>=0.75,1);
select_index=6;
[S,f]=mtspectrumc(analogData_accumulator{select_index,1}(:,s1_index_stim:s2_index_stim)',params);
[S0,f0]=mtspectrumc(analogData0_accumulator{select_index,1}(:,s1_index_stim:s2_index_stim)',params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_2(1,1))
plot(f,(log10(S)),'LineWidth',2,'Color','k');
hold on;
plot(f0,log10(S0),'LineWidth',2,'Color',[0.23,0.51,0.00]);
% ylabel('Power (log(\muV^{2}/Hz))');
ylabel("Raw Power");
% remove Xticks
set(gca,'XTick',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_2(2,1));
s1_index_bl=find(timeVals>=-0.5,1);
s2_index_bl=find(timeVals>=0,1);
[S_bl,f_bl]=mtspectrumc(analogData_accumulator{select_index,1}(:,s1_index_bl:s2_index_bl)',params);
[S0_bl,f0_bl]=mtspectrumc(analogData0_accumulator{select_index,1}(:,s1_index_bl:s2_index_bl)',params);
% plot(f,10*log10(S./S_bl),'LineWidth',2,'Color','k');
plot(f,10*log10(S./S0),'LineWidth',2,'Color','k');
hold on;
% plot(f0,10*log10(S0./S0_bl),'LineWidth',2,'Color',[0.23,0.51,0.00]);
plot(f0,10*log10(S0./S0),'LineWidth',2,'Color',[0.23,0.51,0.00]);
ylabel('Power (dB)');
xlabel('Frequency (Hz)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_3(1,1));
% [h_length_all,p_length_all]=kstest2(length_accumulator_MP{6,1},length_accumulator_hilbert{6,1});
bin_num=16;
min_vals=0;
max_vals=0.8;
edges=linspace(min_vals,max_vals,bin_num+1);
[count,edges]=histcounts(length_accumulator_MP{6,1},edges,'Normalization','probability');
[count2,edges2]=histcounts(length_accumulator_hilbert{6,1},edges,'Normalization','probability');
[count3,edges3]=histcounts(length_accumulator_omp_gear{6,1},edges,'Normalization','probability');
bin_centers_MP=(edges(1:end-1)+edges(2:end))/2;  
bin_centers_hilbert=(edges2(1:end-1)+edges2(2:end))/2;
bin_centers_omp_gear=(edges3(1:end-1)+edges3(2:end))/2;
hold on;
plot(bin_centers_omp_gear, count3, '-', 'LineWidth', 2, 'Color', 'b', 'Marker', '^');
plot(bin_centers_hilbert, count2, '-', 'LineWidth', 2, 'Color', [0.58, 0.0, 0.83], 'Marker', 's');
plot(bin_centers_MP, count, '-', 'LineWidth', 2, 'Color', 'r', 'Marker', 'o');
% plot a vertical dotted line at 0.3 seconds
plot([0.3 0.3], [0 0.8], '--k', 'LineWidth', 1);
xlabel('Burst Length (s)')
ylabel('Fraction of bursts')
legend('OMP-GEAR','Hilbert','MP');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except plotHandles_3 plotHandles_2 plotHandles_1 f
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
subplot(plotHandles_2(1,1))
plot(f_stim,log10(S_stim),'LineWidth',2,'Color','m');
hold on;
% % draw dotted orange lines at 20 and 32 Hz and dotted blue lines at 36 and 65 Hz
color_orange=[0.9, 0.4, 0.0];
% % lighter shade of orange
% % color_orange=[1 0.5 0];
color_blue=[0 0 1];
line([slow_gamma_freq(1) slow_gamma_freq(1)], [min(log10(S_stim))-0.5 max(log10(S_stim))+0.5], 'Color', color_blue, 'LineStyle', '--','LineWidth',2);
line([slow_gamma_freq(2) slow_gamma_freq(2)], [min(log10(S_stim))-0.5 max(log10(S_stim))+0.5], 'Color', color_blue, 'LineStyle', '--','LineWidth',2);
line([fast_gamma_freq(1) fast_gamma_freq(1)], [min(log10(S_stim))-0.5 max(log10(S_stim))+0.5], 'Color', color_orange, 'LineStyle', '--','LineWidth',2);
line([fast_gamma_freq(2) fast_gamma_freq(2)], [min(log10(S_stim))-0.5 max(log10(S_stim))+0.5], 'Color', color_orange, 'LineStyle', '--','LineWidth',2);
ylim([min(log10(S_stim))-0.4 max(log10(S_stim))+0.4]);
legend('Synthetic','Spontaneous','Stimulus-driven');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_2(2,1));
plot(f_stim,10*log10(S_stim./S_bl),'LineWidth',2,'Color','m');
hold on;
line([slow_gamma_freq(1) slow_gamma_freq(1)], [min(10*log10(S_stim./S_bl))-3 max(10*log10(S_stim./S_bl))+3], 'Color', color_blue, 'LineStyle', '--','LineWidth',2);
line([slow_gamma_freq(2) slow_gamma_freq(2)], [min(10*log10(S_stim./S_bl))-3 max(10*log10(S_stim./S_bl))+3], 'Color', color_blue, 'LineStyle', '--','LineWidth',2);
line([fast_gamma_freq(1) fast_gamma_freq(1)], [min(10*log10(S_stim./S_bl))-3 max(10*log10(S_stim./S_bl))+3], 'Color', color_orange, 'LineStyle', '--','LineWidth',2);
line([fast_gamma_freq(2) fast_gamma_freq(2)], [min(10*log10(S_stim./S_bl))-3 max(10*log10(S_stim./S_bl))+3], 'Color', color_orange, 'LineStyle', '--','LineWidth',2);
ylim([min(10*log10(S_stim./S_bl))-2.8 max(10*log10(S_stim./S_bl))+2.8]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set_axis_ticks_fontsize(plotHandles_1,20,16,1);
set_axis_ticks_fontsize(plotHandles_2,20,16,1);
set_axis_ticks_fontsize(plotHandles_2,20,16,2);
set_axis_ticks_fontsize(plotHandles_3,20,16,1);
set_axis_ticks_fontsize(plotHandles_3,20,16,2);
set_axis_ticks_fontsize(plotHandles_3,20,16,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labels = {'A','B'};
x_positions = [0.012];
y_positions = [0.9, 0.52];  % top and bottom rows
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
labels = {'C','D','E'};
x_positions = [0.43];
y_positions = [0.9, 0.6,0.28];  % top and bottom rows
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
    x=x'; y=y';
    SSresid = sum((y - x).^2);        % residuals from y=x
    SStotal = sum((y - mean(y)).^2);  % total variance of y
    R_squared_val = 1 - SSresid / SStotal;
    % RMSE value
    RMSE_val = sqrt(mean((y - x).^2));
end