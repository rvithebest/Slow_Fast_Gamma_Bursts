clc;clear
f=figure;
f.WindowState="Maximized";
plotHandles=getPlotHandles(2,3,[0.08 0.08 0.9 0.9],0.07,0.06,0);
for Monkey_num=1:2
    clearvars -except Monkey_num f plotHandles
    parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
    displayFlag=0;
    stimulusPeriodS=[0.25 0.75];
    baselinePeriodS=[-0.5 0];
    thresholdFraction=0.25;
    num_iterations=120;
    dict_size=2500000;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Monkey_num==1
        % Monkey- alpaH
        load('gamma_duration_alpaH_MP.mat')
        load((fullfile(parent_file_path,'alpaH_info','parameterCombinations.mat')))
        load(fullfile(parent_file_path,'alpaH_info','badTrials.mat'));
        load(fullfile(parent_file_path,'alpaH_info','alpaHMicroelectrodeRFData.mat'));
        LFP_data_file=dir(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H'));
        LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'}));
        LFP_data_file = natsortfiles({LFP_data_file.name});
        slow_gamma_freq=[20 32];
        fast_gamma_freq=[36 65];
        fast_gamma_freq_cont=[20 36];
        gabor_accumulator=gaborInfo_accumulator_alpaH;
        header_accumulator=header_accumulator_alpaH;
    elseif Monkey_num==2
        % Monkey- kesariH
        load('gamma_duration_kesariH_MP.mat')
        load((fullfile(parent_file_path,'kesariH_info','parameterCombinations.mat')))
        load(fullfile(parent_file_path,'kesariH_info','badTrials_kesari.mat'));
        load(fullfile(parent_file_path,'kesariH_info','kesariHMicroelectrodeRFData_Two.mat'));
        LFP_data_file=dir(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H'));
        LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'}));
        LFP_data_file = natsortfiles({LFP_data_file.name});
        slow_gamma_freq=[20 38];
        fast_gamma_freq=[42 65];
        fast_gamma_freq_cont=[20 42];
        gabor_accumulator=gaborInfo_accumulator_kesariH;
        header_accumulator=header_accumulator_kesariH;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('timeVals.mat')
    num_elec=length(current_electrode);
    counter=1;
    length_gatherer_sg=cell(1,num_elec);
    length_gatherer_fg=cell(1,num_elec);
    length_gatherer_bg=cell(1,num_elec);
    onset_gatherer_sg=cell(1,num_elec);
    onset_gatherer_fg=cell(1,num_elec);
    time_gatherer_sg=cell(1,num_elec);
    time_gatherer_fg=cell(1,num_elec);
    freq_gatherer_sg=cell(1,num_elec);
    freq_gatherer_fg=cell(1,num_elec);
    power_gatherer_sg=zeros(1,num_elec);
    power_gatherer_fg=zeros(1,num_elec);
    freq_gatherer_bg=cell(1,num_elec);
    for i=current_electrode
        if Monkey_num==1
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H',LFP_data_file{i}))
        end
        if Monkey_num==2
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H',LFP_data_file{i}))
        end
        % ORI- 157.5 deg and SF- 1cpd : for analysis
        ORI_num=8; SF_num=2;
        trial_temp=parameterCombinations{:,:,:,SF_num,ORI_num};
        trial_temp=setdiff(trial_temp,badTrials);
        data_temp=analogDataDecimated(trial_temp,:);
        gabor_temp=gabor_accumulator{SF_num,ORI_num,counter};
        header_temp=header_accumulator{SF_num,ORI_num,counter};
        %%%%%% Slow gamma burst computation %%%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,slow_gamma_freq);
        power_gatherer_sg(counter)=diffPower;
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        [length_temp_sg,freq_temp_sg,time_center_temp_sg,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        length_temp_all_trials=[];
        onset_temp_all_trials=[];
        freq_temp_all_trials=[];
        time_temp_all_trials=[];
        for ii=1:length(length_temp_sg)
            if isempty(length_temp_sg{ii}')
                continue;
            end
            reject_idx=find((length_temp_sg{ii}')>0.8);
            length_temp_sg{ii}(reject_idx)=[];
            time_center_temp_sg{ii}(reject_idx)=[];
            freq_temp_sg{ii}(reject_idx)=[];
            length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}'];
            freq_temp_all_trials=[freq_temp_all_trials,freq_temp_sg{ii}'];
            time_temp_all_trials=[time_temp_all_trials,time_center_temp_sg{ii}'];
            % Find the first burst (in terms of time) in each trial
            % Onset time (earliest)
            onset_time_temp_trial_sg=time_center_temp_sg{ii}'-((length_temp_sg{ii}')*0.5);
            % onset_time_temp_trial_sg(onset_time_temp_trial_sg<0)=0; % Ensuring onset time is not before stimulus presentation
            onset_idx=(find((onset_time_temp_trial_sg)==min(onset_time_temp_trial_sg)));
            onset_temp_all_trials=[onset_temp_all_trials,onset_time_temp_trial_sg(onset_idx)]; 
        end
        length_gatherer_sg{counter}=length_temp_all_trials;
        freq_gatherer_sg{counter}=freq_temp_all_trials;
        onset_gatherer_sg{counter}=onset_temp_all_trials;
        time_gatherer_sg{counter}=time_temp_all_trials;
        %%%%%% Fast gamma burst computation %%%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
        power_gatherer_fg(counter)=diffPower;
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        [length_temp_fg,freq_temp_fg,time_center_temp_fg,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        length_temp_all_trials=[];
        onset_temp_all_trials=[];
        freq_temp_all_trials=[];
        time_temp_all_trials=[];
        for ii=1:length(length_temp_fg)
            if isempty(length_temp_fg{ii})
                continue;
            end
            reject_idx=find((length_temp_fg{ii}')>0.8);
            length_temp_fg{ii}(reject_idx)=[];
            time_center_temp_fg{ii}(reject_idx)=[];
            freq_temp_fg{ii}(reject_idx)=[];
            length_temp_all_trials=[length_temp_all_trials,length_temp_fg{ii}'];
            freq_temp_all_trials=[freq_temp_all_trials,freq_temp_fg{ii}'];
            time_temp_all_trials=[time_temp_all_trials,time_center_temp_fg{ii}'];
            % Find the first burst (in terms of time) in each trial
            % Onset time (earliest)
            onset_time_temp_trial_fg=time_center_temp_fg{ii}'-((length_temp_fg{ii}')*0.5);
            % onset_time_temp_trial_fg(onset_time_temp_trial_fg<0)=0; % Ensuring onset time is not before stimulus presentation
            onset_idx=(find((onset_time_temp_trial_fg)==min(onset_time_temp_trial_fg)));
            onset_temp_all_trials=[onset_temp_all_trials,onset_time_temp_trial_fg(onset_idx)]; 
        end
        length_gatherer_fg{counter}=length_temp_all_trials;
        freq_gatherer_fg{counter}=freq_temp_all_trials;
        onset_gatherer_fg{counter}=onset_temp_all_trials;
        time_gatherer_fg{counter}=time_temp_all_trials;
        %%%%%%%%%%%%% broad gamma- freq wise computation %%%%%%%%%%%%%%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq_cont);
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        [length_temp_bg,freq_temp_bg,~,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq_cont,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        % lengths from fast gamma are combined with fast_gamma_freq_cont to have broad gamma
        % Note: length_temp_all_trials and freq_temp_all_trials are not reinitialized
        % FG data is concatenated with FG_cont for frequency-wise plot
        for ii=1:length(length_temp_bg)
            if isempty(length_temp_bg{ii})
                continue;
            end
            reject_idx=find((length_temp_bg{ii}')>0.8);
            length_temp_bg{ii}(reject_idx)=[];
            freq_temp_bg{ii}(reject_idx)=[];
            length_temp_all_trials=[length_temp_all_trials,length_temp_bg{ii}'];
            freq_temp_all_trials=[freq_temp_all_trials,freq_temp_bg{ii}'];
        end
        length_gatherer_bg{counter}=length_temp_all_trials;
        freq_gatherer_bg{counter}=freq_temp_all_trials;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        counter=counter+1;
    end
    length_all_elec_bg=[];
    freq_all_elec_bg=[];
    cycles_all_elec_bg=[];
    cycles_all_elec_sg=[];
    cycles_all_elec_fg=[];
    cycles_gatherer_sg=zeros(1,num_elec);
    cycles_gatherer_fg=zeros(1,num_elec);
    median_gatherer_sg=zeros(1,num_elec);
    median_gatherer_fg=zeros(1,num_elec);
    mean_onset_gatherer_sg=zeros(1,num_elec);
    mean_onset_gatherer_fg=zeros(1,num_elec);
    SEM_onset_gatherer_sg=zeros(1,num_elec);
    SEM_onset_gatherer_fg=zeros(1,num_elec);
    mean_time_gatherer_sg=zeros(1,num_elec);
    mean_time_gatherer_fg=zeros(1,num_elec);
    for i=1:num_elec
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        median_gatherer_sg(i)=median(length_gatherer_sg{i});
        median_gatherer_fg(i)=median(length_gatherer_fg{i});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mean_onset_gatherer_sg(i)=mean(onset_gatherer_sg{i});
        mean_onset_gatherer_fg(i)=mean(onset_gatherer_fg{i});
        SEM_onset_gatherer_sg(i)=std(onset_gatherer_sg{i})/sqrt(length(onset_gatherer_sg{i}));
        SEM_onset_gatherer_fg(i)=std(onset_gatherer_fg{i})/sqrt(length(onset_gatherer_fg{i}));
        mean_time_gatherer_sg(i)=median(time_gatherer_sg{i});
        mean_time_gatherer_fg(i)=median(time_gatherer_fg{i});
        % Cycles= Duration of the burst in seconds * Frequency of the burst in Hz
        cycles_temp_sg=length_gatherer_sg{i}.*freq_gatherer_sg{i};
        cycles_temp_fg=length_gatherer_fg{i}.*freq_gatherer_fg{i};
        cycles_gatherer_sg(i)=mean(cycles_temp_sg);
        cycles_gatherer_fg(i)=mean(cycles_temp_fg);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optimized bin number: 7 for both monkeys
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if less than 20 bursts are detected for either slow or fast gamma, then exclude that electrode
        if (length(length_gatherer_sg{i})<20) || (length(length_gatherer_fg{i})<20)
            median_gatherer_sg(i)=0;
            median_gatherer_fg(i)=0;
            power_gatherer_sg(i)=0;
            power_gatherer_fg(i)=0;
            mean_onset_gatherer_sg(i)=0;
            mean_onset_gatherer_fg(i)=0;
            cycles_gatherer_sg(i)=0;
            cycles_gatherer_fg(i)=0;
            mean_time_gatherer_sg(i)=0;
            mean_time_gatherer_fg(i)=0;
            SEM_onset_gatherer_sg(i)=0;
            SEM_onset_gatherer_fg(i)=0;
            continue;
        end
        freq_all_elec_bg=[freq_all_elec_bg,freq_gatherer_bg{i}];
        length_all_elec_bg=[length_all_elec_bg,length_gatherer_bg{i}];
        cycles_all_elec_sg=[cycles_all_elec_sg,cycles_temp_sg];
        cycles_all_elec_fg=[cycles_all_elec_fg,cycles_temp_fg];
        % Cycles= Duration of the burst in seconds * Frequency of the burst in Hz
        cycles_all_elec_bg=[cycles_all_elec_bg,length_gatherer_bg{i}.*freq_gatherer_bg{i}];
    end
    power_gatherer_sg=power_gatherer_sg(power_gatherer_sg~=0);
    power_gatherer_fg=power_gatherer_fg(power_gatherer_fg~=0);
    median_gatherer_sg=median_gatherer_sg(median_gatherer_sg~=0);
    median_gatherer_fg=median_gatherer_fg(median_gatherer_fg~=0);
    mean_onset_gatherer_sg=mean_onset_gatherer_sg(mean_onset_gatherer_sg~=0);
    mean_onset_gatherer_fg=mean_onset_gatherer_fg(mean_onset_gatherer_fg~=0);
    cycles_gatherer_sg=cycles_gatherer_sg(cycles_gatherer_sg~=0);
    cycles_gatherer_fg=cycles_gatherer_fg(cycles_gatherer_fg~=0);
    mean_time_gatherer_sg=mean_time_gatherer_sg(mean_time_gatherer_sg~=0);
    mean_time_gatherer_fg=mean_time_gatherer_fg(mean_time_gatherer_fg~=0);
    SEM_onset_gatherer_fg=SEM_onset_gatherer_fg(SEM_onset_gatherer_fg~=0);
    SEM_onset_gatherer_sg=SEM_onset_gatherer_sg(SEM_onset_gatherer_sg~=0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    power_gatherer_sg=10*log10(power_gatherer_sg);
    power_gatherer_fg=10*log10(power_gatherer_fg);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,1))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on;
    comparator=0;
    num_elec_2=length(mean_onset_gatherer_sg);
    orange_color=[1 0.5 0];
    for i=1:num_elec_2
       patch_width=0.7;
       x1= i - patch_width/2;
       x2= i + patch_width/2;
       y1_sg= mean_onset_gatherer_sg(i) - SEM_onset_gatherer_sg(i);
       y2_sg= mean_onset_gatherer_sg(i) + SEM_onset_gatherer_sg(i);
       y1_fg= mean_onset_gatherer_fg(i) - SEM_onset_gatherer_fg(i);
       y2_fg= mean_onset_gatherer_fg(i) + SEM_onset_gatherer_fg(i);
       plot(i,mean_onset_gatherer_sg(i),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','b','Color',[0 0 1],'MarkerSize',5); hold on;
       plot(i,mean_onset_gatherer_fg(i),'o','MarkerFaceColor',orange_color,'MarkerEdgeColor',orange_color,'Color',orange_color,'MarkerSize',5); hold on;
       patch([x1 x1 x2 x2], [y1_sg y2_sg y2_sg y1_sg], [0 0 1],'EdgeColor','none','FaceAlpha',0.2); hold on;
       patch([x1 x1 x2 x2], [y1_fg y2_fg y2_fg y1_fg], orange_color,'EdgeColor','none','FaceAlpha',0.2); hold on;
       if mean_onset_gatherer_sg(i) > mean_onset_gatherer_fg(i)
           comparator = comparator + 1;
       end
       ylim([0.1 0.4]);
   end
   legend('Slow \gamma', 'Fast \gamma');
   display(['Slow gamma bursts are delayed than fast gamma bursts in ', num2str(comparator), ' out of ', num2str(num_elec_2), ' electrodes']);
   xlabel('Electrode number');
   ylabel('Onset time (s)');
   % %%%%%%%% Horizontal plot %%%%%%%%%%%%%%%%%%%%%%%%%%%
   view([90 -90]); % Rotates the axes
   set(gca, 'YDir', 'normal'); % Ensure correct label orientation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % violin_swarm_plot_paired(mean_time_gatherer_sg,mean_time_gatherer_fg,1);
    % violin_swarm_plot(cycles_gatherer_sg,cycles_gatherer_fg);
    % ylabel("Median No. of cycles");
    % Earlier option: histogram of time centers
    % histogram_all_burst(all_time_gatherer_sg,all_time_gatherer_fg,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,2))
    violin_swarm_plot_paired(mean_onset_gatherer_sg,mean_onset_gatherer_fg,2);
    % ylim([0.25 0.75]);
    ylim([0.1 0.4]);
    ylabel("Mean Onset Time");
    if Monkey_num==1
        save('Onset_center_MP_fig_4_M1.mat','mean_onset_gatherer_sg','mean_onset_gatherer_fg');
    elseif Monkey_num==2
        save('Onset_center_MP_fig_4_M2.mat','mean_onset_gatherer_sg','mean_onset_gatherer_fg');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,3))
    hold on;
    bin_size=4; % from 20 to 65 Hz
    % find corresponding bin for each frequency and thus corresponding length
    freq_bins=20:bin_size:65;
    freq_bin_centers=freq_bins(1:end-1)+bin_size/2;
    freq_bins_num=discretize(freq_all_elec_bg,freq_bins);
    % find median length and SEM corresponding to each bin
    median_length_freq_bins=zeros(1,length(freq_bin_centers));
    sem_length_freq_bins=zeros(1,length(freq_bin_centers));
    median_cycles_freq_bins=zeros(1,length(freq_bin_centers));
    sem_cycles_freq_bins=zeros(1,length(freq_bin_centers));
    for i=1:(length(freq_bins)-1)
        temp_indices=find(freq_bins_num==i);
        temp_lengths=length_all_elec_bg(temp_indices);
        temp_cycles=cycles_all_elec_bg(temp_indices);
        median_length_freq_bins(i)=median(temp_lengths);
        sem_length_freq_bins(i)=getSEMedian(temp_lengths,1000);
        median_cycles_freq_bins(i)=median(temp_cycles);
        sem_cycles_freq_bins(i)=getSEMedian(temp_cycles,1000);
    end
    errorbar(freq_bin_centers,median_length_freq_bins,sem_length_freq_bins,'-o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',1.5);
    if Monkey_num==1
        % Label the first three frequency centers data points with blue color (slow gamma)
        plot(freq_bin_centers(1:3), median_length_freq_bins(1:3), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1], 'Color', [0 0 1], 'LineWidth', 1.5);
        plot(freq_bin_centers(5:11), median_length_freq_bins(5:11), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0], 'Color', [1 0.5 0], 'LineWidth', 1.5);
    elseif Monkey_num==2
        plot(freq_bin_centers(1:4), median_length_freq_bins(1:4), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1], 'Color', [0 0 1], 'LineWidth', 1.5);
        plot(freq_bin_centers(6:11), median_length_freq_bins(6:11), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0], 'Color', [1 0.5 0], 'LineWidth', 1.5);
    end
    xlabel('Frequency (Hz)');
    ylabel('Burst length (s)');
    set_axis_ticks_fontsize(plotHandles,22,16,Monkey_num);
    %%%% Cycles plot %%%%%%
    % figure;
    % errorbar(freq_bin_centers,median_cycles_freq_bins,sem_cycles_freq_bins,'-o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',1.5);
end
labels = {'A','B','C','D','E','F'};
x_positions = [0.03, 0.36 0.68];
y_positions = [0.88, 0.41];  % top and bottom rows
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annotation('textbox',...
[0.05 0.68 0.08 0.09],...
'String',{'Monkey 1'},...
'Rotation',90,...
'FontWeight','bold',...
'FontSize',20,...
'FontName','Helvetica',...
'EdgeColor',[1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 annotation('textbox',...
[0.05 0.21 0.08 0.09],...
'String',{'Monkey 2'},...
'Rotation',90,...
'FontWeight','bold',...
'FontSize',20,...
'FontName','Helvetica',...
'EdgeColor',[1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annotation('textbox', ...
    [0.06, 0.87, 0.1, 0.1], ...
    'String', 'N=77', ...
    'FontSize', 28, ...
    'FontWeight', 'Bold', ...
    'EdgeColor', 'none', ...
    'FontName', 'Helvetica');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annotation('textbox', ...
    [0.06, 0.4, 0.1, 0.1], ...
    'String', 'N=31', ...
    'FontSize', 28, ...
    'FontWeight', 'Bold', ...
    'EdgeColor', 'none', ...
    'FontName', 'Helvetica');