clc;clear;close all
f=figure;
f.WindowState="Maximized";
plotHandles=getPlotHandles(2,4,[0.08 0.08 0.9 0.85],0.07,0.06,0);
for Monkey_num=1:2
    clearvars -except Monkey_num f plotHandles
    parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
    displayFlag=0;
    stimulusPeriodS=[0.25 0.75];
    baselinePeriodS=[-0.5 0];
    thresholdFraction=0.5;
    filterOrder=4;
    num_iterations=120;
    dict_size=2500000;
    if Monkey_num==1
      % Monkey- alpaH
        load('gamma_duration_alpaH_OMP_GEAR.mat');
        load((fullfile(parent_file_path,'alpaH_info','parameterCombinations.mat')))
        load(fullfile(parent_file_path,'alpaH_info','badTrials.mat'));
        load(fullfile(parent_file_path,'alpaH_info','alpaHMicroelectrodeRFData.mat'));
        LFP_data_file=dir(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H'));
        LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'}));
        LFP_data_file = natsortfiles({LFP_data_file.name});
        slow_gamma_freq=[20 32];
        fast_gamma_freq=[36 65];
        gabor_accumulator=gaborInfo_accumulator_ORI;
        header_accumulator=header_accumulator_ORI;
    elseif Monkey_num==2
        % Monkey- kesariH
        load('gamma_duration_kesariH_OMP_GEAR.mat');
        load((fullfile(parent_file_path,'kesariH_info','parameterCombinations.mat')))
        load(fullfile(parent_file_path,'kesariH_info','badTrials_kesari.mat'));
        load(fullfile(parent_file_path,'kesariH_info','kesariHMicroelectrodeRFData_Two.mat'));
        LFP_data_file=dir(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H'));
        LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'}));
        LFP_data_file = natsortfiles({LFP_data_file.name});
        slow_gamma_freq=[20 38];
        fast_gamma_freq=[42 65];
        gabor_accumulator=gaborInfo_accumulator_ORI;
        header_accumulator=header_accumulator_ORI;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('timeVals.mat')
    num_elec=length(current_electrode);
    counter=1;
    length_gatherer_sg_hilbert=cell(1,num_elec);
    length_gatherer_fg_hilbert=cell(1,num_elec);
    length_gatherer_sg_OMP_gear=cell(1,num_elec);
    length_gatherer_fg_OMP_gear=cell(1,num_elec);
    %%% for power gatherer
    power_gatherer_sg_hilbert=zeros(1,num_elec);
    power_gatherer_fg_hilbert=zeros(1,num_elec);
    power_gatherer_sg_OMP_gear=zeros(1,num_elec);
    power_gatherer_fg_OMP_gear=zeros(1,num_elec);
    % for onset gatherer
    onset_gatherer_sg_hilbert=cell(1,num_elec);
    onset_gatherer_fg_hilbert=cell(1,num_elec);
    onset_gatherer_sg_OMP_gear=cell(1,num_elec);
    onset_gatherer_fg_OMP_gear=cell(1,num_elec);
    %%% for se start karna hai
    for i=current_electrode
        if Monkey_num==1
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H',LFP_data_file{i}))
            gabor_temp=gabor_accumulator{1,1,counter};
            header_temp=header_accumulator{1,1,counter};
        end
        if Monkey_num==2
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H',LFP_data_file{i}))
            gabor_temp=gabor_accumulator{1,8,counter};
            header_temp=header_accumulator{1,8,counter};
        end
        % ORI- 157.5 deg and SF- 1cpd : for analysis
        ORI_num=8; SF_num=2;
        trial_temp=parameterCombinations{:,:,:,SF_num,ORI_num};
        trial_temp=setdiff(trial_temp,badTrials);
        data_temp=analogDataDecimated(trial_temp,:);
        %%%%%% Slow gamma burst computation %%%%%%
        %%% Hilbert transform %%%%%%%%%%%%%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,slow_gamma_freq);
        power_gatherer_sg_hilbert(counter)=diffPower;
        power_gatherer_sg_OMP_gear(counter)=diffPower;
        thresholdFactor=(thresholdFraction*diffPower);
        [length_temp_sg,time_center_temp_sg]= getBurstLengthHilbert(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,filterOrder);
        length_temp_all_trials=[];
        onset_temp_all_trials=[];
        for ii=1:length(length_temp_sg)
            if isempty(length_temp_sg{ii})
                continue;
            end
            reject_idx=find((length_temp_sg{ii})>0.8);
            length_temp_sg{ii}(reject_idx)=[];
            time_center_temp_sg{ii}(reject_idx)=[];
            length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}];
            % Find the first burst (in terms of time) in each trial
            % Onset time (earliest)
            onset_time_temp_trial_sg=time_center_temp_sg{ii}-((length_temp_sg{ii})*0.5);
            % onset_time_temp_trial_sg(onset_time_temp_trial_sg<0)=0; % Ensuring onset time is not before stimulus presentation
            onset_idx=(find((onset_time_temp_trial_sg)==min(onset_time_temp_trial_sg)));
            onset_temp_all_trials=[onset_temp_all_trials,onset_time_temp_trial_sg(onset_idx)]; 
        end
        length_gatherer_sg_hilbert{counter}=length_temp_all_trials;
        onset_gatherer_sg_hilbert{counter}=onset_temp_all_trials;
        %%%%%% OMP-GEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thresholdFactor=sqrt(0.5*diffPower);
        [length_temp_sg,~,time_center_temp_sg,~,~,~]= getBurstLength_all(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp,'OMP-GEAR');
        length_temp_all_trials=[];
        onset_temp_all_trials=[];
        for ii=1:length(length_temp_sg)
            if isempty(length_temp_sg{ii}')
                continue;
            end
            reject_idx=find((length_temp_sg{ii}')>0.8);
            length_temp_sg{ii}(reject_idx)=[];
            time_center_temp_sg{ii}(reject_idx)=[];
            length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}'];
            % Find the first burst (in terms of time) in each trial
            % Onset time (earliest)
            onset_time_temp_trial_sg=(time_center_temp_sg{ii}')-((length_temp_sg{ii}')*0.5);
            % onset_time_temp_trial_sg(onset_time_temp_trial_sg<0)=0; % Ensuring onset time is not before stimulus presentation
            onset_idx=(find((onset_time_temp_trial_sg)==min(onset_time_temp_trial_sg)));
            onset_temp_all_trials=[onset_temp_all_trials,onset_time_temp_trial_sg(onset_idx)];
        end
        length_gatherer_sg_OMP_gear{counter}=length_temp_all_trials;
        onset_gatherer_sg_OMP_gear{counter}=onset_temp_all_trials;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Fast gamma burst computation %%%%%%
        %%% Hilbert transform %%%%%%%%%%%%%%%%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
        power_gatherer_fg_hilbert(counter)=diffPower;
        power_gatherer_fg_OMP_gear(counter)=diffPower;
        power_gatherer_fg_feingold(counter)=diffPower;
        power_gatherer_fg_wavelet(counter)=diffPower;
        thresholdFactor=(thresholdFraction*diffPower);
        [length_temp_fg,time_center_temp_fg]= getBurstLengthHilbert(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,filterOrder);
        length_temp_all_trials=[];
        onset_temp_all_trials=[];
        for ii=1:length(length_temp_fg)
            if isempty(length_temp_fg{ii})
                continue;
            end
            reject_idx=find((length_temp_fg{ii})>0.8);
            length_temp_fg{ii}(reject_idx)=[];
            time_center_temp_fg{ii}(reject_idx)=[];
            length_temp_all_trials=[length_temp_all_trials,length_temp_fg{ii}];
            % Find the first burst (in terms of time) in each trial
            % Onset time (earliest)
            onset_time_temp_trial_fg=(time_center_temp_fg{ii}')-((length_temp_fg{ii}')*0.5);
            % onset_time_temp_trial_fg(onset_time_temp_trial_fg<0)=0; % Ensuring onset time is not before stimulus presentation
            onset_idx=(find((onset_time_temp_trial_fg)==min(onset_time_temp_trial_fg)));
            onset_temp_all_trials=[onset_temp_all_trials,onset_time_temp_trial_fg(onset_idx)];
        end
        length_gatherer_fg_hilbert{counter}=length_temp_all_trials;
        onset_gatherer_fg_hilbert{counter}=onset_temp_all_trials;
        %%%%%% OMP-GEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thresholdFactor=sqrt(0.5*diffPower);
        [length_temp_fg,~,time_center_temp_fg,~,~,~]= getBurstLength_all(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp,'OMP-GEAR');
        length_temp_all_trials=[];
        onset_temp_all_trials=[];
        for ii=1:length(length_temp_fg)
            if isempty(length_temp_fg{ii})
                continue;
            end
            reject_idx=find((length_temp_fg{ii})>0.8);
            length_temp_fg{ii}(reject_idx)=[];
            time_center_temp_fg{ii}(reject_idx)=[];
            length_temp_all_trials=[length_temp_all_trials,length_temp_fg{ii}'];
            % Find the first burst (in terms of time) in each trial
            % Onset time (earliest)
            onset_time_temp_trial_fg=(time_center_temp_fg{ii}')-((length_temp_fg{ii}')*0.5);
            % onset_time_temp_trial_fg(onset_time_temp_trial_fg<0)=0; % Ensuring onset time is not before stimulus presentation
            onset_idx=(find((onset_time_temp_trial_fg)==min(onset_time_temp_trial_fg)));
            onset_temp_all_trials=[onset_temp_all_trials,onset_time_temp_trial_fg(onset_idx)];
        end
        length_gatherer_fg_OMP_gear{counter}=length_temp_all_trials;
        onset_gatherer_fg_OMP_gear{counter}=onset_temp_all_trials;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        counter=counter+1;
    end
    length_all_elec_sg_hilbert=[];
    length_all_elec_fg_hilbert=[];
    median_gatherer_sg_hilbert=zeros(1,num_elec);
    median_gatherer_fg_hilbert=zeros(1,num_elec);
    length_all_elec_sg_OMP_gear=[];
    length_all_elec_fg_OMP_gear=[];
    median_gatherer_sg_OMP_gear=zeros(1,num_elec);
    median_gatherer_fg_OMP_gear=zeros(1,num_elec);
    mean_onset_gatherer_sg_OMP_gear=zeros(1,num_elec);
    mean_onset_gatherer_fg_OMP_gear=zeros(1,num_elec);
    mean_onset_gatherer_sg_hilbert=zeros(1,num_elec);
    mean_onset_gatherer_fg_hilbert=zeros(1,num_elec);
    for i=1:num_elec
        %%%%%%%%%%%%% Hilbert transform %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        median_gatherer_sg_hilbert(i)=median(length_gatherer_sg_hilbert{i});
        median_gatherer_fg_hilbert(i)=median(length_gatherer_fg_hilbert{i});
        mean_onset_gatherer_sg_hilbert(i)=mean(onset_gatherer_sg_hilbert{i});
        mean_onset_gatherer_fg_hilbert(i)=mean(onset_gatherer_fg_hilbert{i});
        if (length(length_gatherer_sg_hilbert{i})<20 || length(length_gatherer_fg_hilbert{i})<20)
            median_gatherer_sg_hilbert(i)=0;
            median_gatherer_fg_hilbert(i)=0;
            power_gatherer_sg_hilbert(i)=0;
            power_gatherer_fg_hilbert(i)=0;
            mean_onset_gatherer_sg_hilbert(i)=0;
            mean_onset_gatherer_fg_hilbert(i)=0;
        end
        %%%%% OMP-GEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        median_gatherer_sg_OMP_gear(i)=median(length_gatherer_sg_OMP_gear{i});
        median_gatherer_fg_OMP_gear(i)=median(length_gatherer_fg_OMP_gear{i});
        mean_onset_gatherer_sg_OMP_gear(i)=mean(onset_gatherer_sg_OMP_gear{i});
        mean_onset_gatherer_fg_OMP_gear(i)=mean(onset_gatherer_fg_OMP_gear{i});
        if (length(length_gatherer_sg_OMP_gear{i})<20) || (length(length_gatherer_fg_OMP_gear{i})<20)
            median_gatherer_sg_OMP_gear(i)=0;
            median_gatherer_fg_OMP_gear(i)=0;
            power_gatherer_sg_OMP_gear(i)=0;
            power_gatherer_fg_OMP_gear(i)=0;
            mean_onset_gatherer_sg_OMP_gear(i)=0;
            mean_onset_gatherer_fg_OMP_gear(i)=0;
        end
    end
    power_gatherer_sg_hilbert=power_gatherer_sg_hilbert(power_gatherer_sg_hilbert~=0);
    power_gatherer_fg_hilbert=power_gatherer_fg_hilbert(power_gatherer_fg_hilbert~=0);
    median_gatherer_sg_hilbert=median_gatherer_sg_hilbert(median_gatherer_sg_hilbert~=0);
    median_gatherer_fg_hilbert=median_gatherer_fg_hilbert(median_gatherer_fg_hilbert~=0);
    power_gatherer_sg_OMP_gear=power_gatherer_sg_OMP_gear(power_gatherer_sg_OMP_gear~=0);
    power_gatherer_fg_OMP_gear=power_gatherer_fg_OMP_gear(power_gatherer_fg_OMP_gear~=0);
    median_gatherer_sg_OMP_gear=median_gatherer_sg_OMP_gear(median_gatherer_sg_OMP_gear~=0);
    median_gatherer_fg_OMP_gear=median_gatherer_fg_OMP_gear(median_gatherer_fg_OMP_gear~=0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    power_gatherer_sg_hilbert=10*log10(power_gatherer_sg_hilbert);
    power_gatherer_fg_hilbert=10*log10(power_gatherer_fg_hilbert);
    power_gatherer_sg_OMP_gear=10*log10(power_gatherer_sg_OMP_gear);
    power_gatherer_fg_OMP_gear=10*log10(power_gatherer_fg_OMP_gear);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mean_onset_gatherer_sg_hilbert=mean_onset_gatherer_sg_hilbert(mean_onset_gatherer_sg_hilbert~=0);
    mean_onset_gatherer_fg_hilbert=mean_onset_gatherer_fg_hilbert(mean_onset_gatherer_fg_hilbert~=0);
    mean_onset_gatherer_sg_OMP_gear=mean_onset_gatherer_sg_OMP_gear(mean_onset_gatherer_sg_OMP_gear~=0);
    mean_onset_gatherer_fg_OMP_gear=mean_onset_gatherer_fg_OMP_gear(mean_onset_gatherer_fg_OMP_gear~=0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Hilbert transform %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,1))
    if Monkey_num==1
       load('length_hilbert_matched_power_M1.mat','matched_sg_lengths','matched_fg_lengths','matched_sg_indices','matched_fg_indices');
    end
    if Monkey_num==2
        load('length_hilbert_matched_power_M2.mat','matched_sg_lengths','matched_fg_lengths','matched_sg_indices','matched_fg_indices');
    end
    % [matched_sg_indices,matched_fg_indices]=power_matching_hist(power_gatherer_sg_hilbert,power_gatherer_fg_hilbert);
    matched_sg_lengths=median_gatherer_sg_hilbert(matched_sg_indices);
    matched_fg_lengths=median_gatherer_fg_hilbert(matched_fg_indices);
    violin_swarm_plot(matched_sg_lengths,matched_fg_lengths);
    if Monkey_num==1
        title('Hilbert','FontSize',24,'FontWeight','bold','FontName','Helvetica');
    end
    subplot(plotHandles(Monkey_num,2))
    violin_swarm_plot_paired(mean_onset_gatherer_sg_hilbert, mean_onset_gatherer_fg_hilbert,2);
    if Monkey_num==1
       save("Onset_gather_hilbert_M1.mat","mean_onset_gatherer_sg_hilbert","mean_onset_gatherer_fg_hilbert");
    elseif Monkey_num==2
       save("Onset_gather_hilbert_M2.mat","mean_onset_gatherer_sg_hilbert","mean_onset_gatherer_fg_hilbert");
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% OMP-GEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,3))
    if Monkey_num==1
       ylim([0 0.6])
       load('length_OMP_gear_matched_power_M1.mat','matched_sg_lengths','matched_fg_lengths','matched_sg_indices','matched_fg_indices');
    end 
    if Monkey_num==2
        ylim([0 0.4])
        load('length_OMP_gear_matched_power_M2.mat','matched_sg_lengths','matched_fg_lengths','matched_sg_indices','matched_fg_indices');
    end
    % [matched_sg_indices,matched_fg_indices]=power_matching_hist(power_gatherer_sg_OMP_gear,power_gatherer_fg_OMP_gear);
    %  matched_sg_lengths=median_gatherer_sg_OMP_gear(matched_sg_indices);
    %  matched_fg_lengths=median_gatherer_fg_OMP_gear(matched_fg_indices);
    violin_swarm_plot(matched_sg_lengths,matched_fg_lengths);
    if Monkey_num==1
        title('OMP-GEAR','FontSize',24,'FontWeight','bold','FontName','Helvetica');
    end
    subplot(plotHandles(Monkey_num,4))
    violin_swarm_plot_paired(mean_onset_gatherer_sg_OMP_gear, mean_onset_gatherer_fg_OMP_gear,2)
    if Monkey_num==1
       save("Onset_gather_OMP_gear_M1.mat","mean_onset_gatherer_sg_OMP_gear","mean_onset_gatherer_fg_OMP_gear");
    elseif Monkey_num==2
       save("Onset_gather_OMP_gear_M2.mat","mean_onset_gatherer_sg_OMP_gear","mean_onset_gatherer_fg_OMP_gear");
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set_axis_ticks_fontsize(plotHandles,22,16,Monkey_num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labels = {'A','B','C','D','E','F','G','H'};
x_positions = [0.03, 0.27, 0.51, 0.76];
y_positions = [0.9, 0.42];  % top and bottom rows
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