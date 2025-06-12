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
    thresholdFraction=0.5;
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
    time_gatherer_sg=cell(1,num_elec);
    time_gatherer_fg=cell(1,num_elec);
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
        [length_temp_sg,~,time_temp_sg,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        length_temp_all_trials=[];
        time_temp_all_trials=[];
        for ii=1:length(length_temp_sg)
            length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}'];
            time_temp_all_trials=[time_temp_all_trials,time_temp_sg{ii}'];
        end
        length_gatherer_sg{counter}=length_temp_all_trials;
        time_gatherer_sg{counter}=time_temp_all_trials;
        %%%%%% Fast gamma burst computation %%%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
        power_gatherer_fg(counter)=diffPower;
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        [length_temp_fg,freq_temp_fg,time_temp_fg,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        length_temp_all_trials=[];
        time_temp_all_trials=[];
        freq_temp_all_trials=[];
        for ii=1:length(length_temp_fg)
            length_temp_all_trials=[length_temp_all_trials,length_temp_fg{ii}'];
            freq_temp_all_trials=[freq_temp_all_trials,freq_temp_fg{ii}'];
            time_temp_all_trials=[time_temp_all_trials,time_temp_fg{ii}'];
        end
        length_gatherer_fg{counter}=length_temp_all_trials;
        time_gatherer_fg{counter}=time_temp_all_trials;
        %%%% broad gamma- freq wise computation %%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq_cont);
        thresholdFactor=thresholdFraction*sqrt(diffPower);
        [length_temp_bg,freq_temp_bg,~,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq_cont,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        % lengths from fast gamma are combined with fast_gamma_freq_cont to have broad gamma
        for ii=1:length(length_temp_bg)
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
    median_gatherer_sg=zeros(1,num_elec);
    median_gatherer_fg=zeros(1,num_elec);
    median_time_gatherer_sg=zeros(1,num_elec);
    median_time_gatherer_fg=zeros(1,num_elec);
    all_time_gatherer_sg=[];
    all_time_gatherer_fg=[];
    for i=1:num_elec
        reject_indices_sg=find(length_gatherer_sg{i}>0.8);
        reject_indices_fg=find(length_gatherer_fg{i}>0.8);
        length_gatherer_sg{i}(reject_indices_sg)=[];
        length_gatherer_fg{i}(reject_indices_fg)=[];
        time_gatherer_sg{i}(reject_indices_sg)=[];
        time_gatherer_fg{i}(reject_indices_fg)=[];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        reject_sg_indices=find(length_gatherer_sg{i} < 0.1);
        reject_fg_indices=find(length_gatherer_fg{i} < 0.1);
        time_gatherer_sg{i}(reject_sg_indices)=[];
        time_gatherer_fg{i}(reject_fg_indices)=[];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        median_gatherer_sg(i)=median(length_gatherer_sg{i});
        median_gatherer_fg(i)=median(length_gatherer_fg{i});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        median_time_gatherer_sg(i)=median(time_gatherer_sg{i});
        median_time_gatherer_fg(i)=median(time_gatherer_fg{i});
        % Plot SEM of median time gatherer
        SEM_time_gatherer_sg(i)=getSEMedian(time_gatherer_sg{i},1000);
        SEM_time_gatherer_fg(i)=getSEMedian(time_gatherer_fg{i},1000);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        all_time_gatherer_sg=[all_time_gatherer_sg,time_gatherer_sg{i}];
        all_time_gatherer_fg=[all_time_gatherer_fg,time_gatherer_fg{i}];
        % Optimized bin number: 7 for both monkeys
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        reject_indices_bg=find(length_gatherer_bg{i}>0.8);
        length_gatherer_bg{i}(reject_indices_bg)=[];
        freq_gatherer_bg{i}(reject_indices_bg)=[];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % if less than  median length is less than 100 msec for either slow or fast gamma, then exclude that electrode
        % if median_gatherer_sg(i)<0.1 || median_gatherer_fg(i)<0.1
        %     median_gatherer_sg(i)=0;
        %     median_gatherer_fg(i)=0;
        %     power_gatherer_sg(i)=0;
        %     power_gatherer_fg(i)=0;
        %     median_time_gatherer_sg(i)=0;
        %     median_time_gatherer_fg(i)=0;
        %     SEM_time_gatherer_sg(i)=0;
        %     SEM_time_gatherer_fg(i)=0;
        %     continue;
        % end
        % if less than 20 bursts are detected for either slow or fast gamma, then exclude that electrode
        if length(length_gatherer_sg{i})<20 || length(length_gatherer_fg{i})<20
            median_gatherer_sg(i)=0;
            median_gatherer_fg(i)=0;
            power_gatherer_sg(i)=0;
            power_gatherer_fg(i)=0;
            median_time_gatherer_sg(i)=0;
            median_time_gatherer_fg(i)=0;
            SEM_time_gatherer_sg(i)=0;
            SEM_time_gatherer_fg(i)=0;
            continue;
        end
        freq_all_elec_bg=[freq_all_elec_bg,freq_gatherer_bg{i}];
        length_all_elec_bg=[length_all_elec_bg,length_gatherer_bg{i}];
    end
    power_gatherer_sg=power_gatherer_sg(power_gatherer_sg~=0);
    power_gatherer_fg=power_gatherer_fg(power_gatherer_fg~=0);
    median_gatherer_sg=median_gatherer_sg(median_gatherer_sg~=0);
    median_gatherer_fg=median_gatherer_fg(median_gatherer_fg~=0);
    median_time_gatherer_sg=median_time_gatherer_sg(median_time_gatherer_sg~=0);
    median_time_gatherer_fg=median_time_gatherer_fg(median_time_gatherer_fg~=0);
    SEM_time_gatherer_sg=SEM_time_gatherer_sg(SEM_time_gatherer_sg~=0);
    SEM_time_gatherer_fg=SEM_time_gatherer_fg(SEM_time_gatherer_fg~=0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    power_gatherer_sg=10*log10(power_gatherer_sg);
    power_gatherer_fg=10*log10(power_gatherer_fg);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,1))
    hold on;
    comparator=0;
    num_elec_2=length(median_time_gatherer_sg);
    orange_color=[1 0.5 0];
    for i = 1:num_elec_2
        %plot((median_time_gatherer_sg(i)), i,'w','Marker','o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','b'); hold on;
        % errorbar(i,median_time_gatherer_sg(i), SEM_time_gatherer_sg(i), 'w', 'Marker', 'o', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 1], 'Color', [0 0 1], 'LineWidth', 1.5); hold on;
        %plot((median_time_gatherer_fg(i)), i,'w','Marker','o','MarkerFaceColor',orange_color,'MarkerEdgeColor',orange_color); hold on;
        % errorbar(i,median_time_gatherer_fg(i), SEM_time_gatherer_fg(i), 'w', 'Marker', 'o', 'MarkerFaceColor', orange_color, 'MarkerEdgeColor', orange_color, 'Color', orange_color, 'LineWidth', 1.5); hold on;
        patch_width=0.5;
        x1= i - patch_width/2;
        x2= i + patch_width/2;
        y1_sg= median_time_gatherer_sg(i) - SEM_time_gatherer_sg(i);
        y2_sg= median_time_gatherer_sg(i) + SEM_time_gatherer_sg(i);
        y1_fg= median_time_gatherer_fg(i) - SEM_time_gatherer_fg(i);
        y2_fg= median_time_gatherer_fg(i) + SEM_time_gatherer_fg(i);
        plot(i,median_time_gatherer_sg(i),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','b','Color',[0 0 1],'MarkerSize',5); hold on;
        plot(i,median_time_gatherer_fg(i),'o','MarkerFaceColor',orange_color,'MarkerEdgeColor',orange_color,'Color',orange_color,'MarkerSize',5); hold on;
        patch([x1 x1 x2 x2], [y1_sg y2_sg y2_sg y1_sg], [0 0 1],'EdgeColor','none','FaceAlpha',0.2); hold on;
        patch([x1 x1 x2 x2], [y1_fg y2_fg y2_fg y1_fg], orange_color,'EdgeColor','none','FaceAlpha',0.2); hold on;
        if median_time_gatherer_sg(i) > median_time_gatherer_fg(i)
            comparator = comparator + 1;
        end
        ylim([0.25 0.75]);
    end
    legend('Slow \gamma', 'Fast \gamma');
    display(['Slow gamma bursts are delayed than fast gamma bursts in ', num2str(comparator), ' out of ', num2str(num_elec_2), ' electrodes']);
    xlabel('Electrode number');
    ylabel('Time center (s)');
   %%%%%%%% Horizontal plot %%%%%%%%%%%%%%%%%%%%%%%%%%%
    view([90 -90]); % Rotates the axes
    set(gca, 'YDir', 'normal'); % Ensure correct label orientation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ylabel('Electrode number');
    % xlabel('Median time center (s)');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,2))
    violin_swarm_plot_paired(median_time_gatherer_sg,median_time_gatherer_fg);
    ylim([0.25 0.75]);
    if Monkey_num==1
        save('time_center_MP_fig_3_M1.mat','median_time_gatherer_sg','median_time_gatherer_fg');
    elseif Monkey_num==2
        save('time_center_MP_fig_3_M2.mat','median_time_gatherer_sg','median_time_gatherer_fg');
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
    for i=1:(length(freq_bins)-1)
        temp_indices=find(freq_bins_num==i);
        temp_lengths=length_all_elec_bg(temp_indices);
        median_length_freq_bins(i)=median(temp_lengths);
        sem_length_freq_bins(i)=getSEMedian(temp_lengths);
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
end