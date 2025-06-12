clc;clear;close all
f=figure;
f.WindowState="Maximized";
plotHandles=getPlotHandles(2,3,[0.08 0.08 0.9 0.85],0.07,0.06,0);
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
    length_gatherer_sg_feingold=cell(1,num_elec);
    length_gatherer_fg_feingold=cell(1,num_elec);
    length_gatherer_sg_wavelet=cell(1,num_elec);
    length_gatherer_fg_wavelet=cell(1,num_elec);
    %%% for power gatherer
    power_gatherer_sg_hilbert=zeros(1,num_elec);
    power_gatherer_fg_hilbert=zeros(1,num_elec);
    power_gatherer_sg_OMP_gear=zeros(1,num_elec);
    power_gatherer_fg_OMP_gear=zeros(1,num_elec);
    power_gatherer_sg_feingold=zeros(1,num_elec);
    power_gatherer_fg_feingold=zeros(1,num_elec);
    power_gatherer_sg_wavelet=zeros(1,num_elec);
    power_gatherer_fg_wavelet=zeros(1,num_elec);
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
        power_gatherer_sg_feingold(counter)=diffPower;
        power_gatherer_sg_wavelet(counter)=diffPower;
        thresholdFactor=(thresholdFraction*diffPower);
        [length_temp_sg,~]= getBurstLengthHilbert(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,filterOrder);
        length_temp_all_trials=[];
        for ii=1:length(length_temp_sg)
             length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}];
        end
        length_gatherer_sg_hilbert{counter}=length_temp_all_trials;
        %%%%%% Wavelet transform %%%%%%%%%%%%%%%%
        [length_temp_sg,~]= getBurstLengthWavelet(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq);
        length_temp_all_trials=[];
        for ii=1:length(length_temp_sg)
            length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}];
        end
        length_gatherer_sg_wavelet{counter}=length_temp_all_trials;
        %%%%%% Feingold %%%%%%%%%%%%%%%%
        [length_temp_sg]= getBurstLengthFeingold(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,filterOrder);
        length_temp_all_trials=[];
        for ii=1:length(length_temp_sg)
            length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}];
        end
        length_gatherer_sg_feingold{counter}=length_temp_all_trials;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% OMP-GEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thresholdFactor=0.75*sqrt(diffPower);
        [length_temp_sg,~,~,~,~,~]= getBurstLength_all(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp,'OMP-GEAR');
        length_temp_all_trials=[];
        for ii=1:length(length_temp_sg)
            length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}'];
        end
        length_gatherer_sg_OMP_gear{counter}=length_temp_all_trials;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Fast gamma burst computation %%%%%%
        %%% Hilbert transform %%%%%%%%%%%%%%%%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
        power_gatherer_fg_hilbert(counter)=diffPower;
        power_gatherer_fg_OMP_gear(counter)=diffPower;
        power_gatherer_fg_feingold(counter)=diffPower;
        power_gatherer_fg_wavelet(counter)=diffPower;
        thresholdFactor=(thresholdFraction*diffPower);
        [length_temp_fg,~]= getBurstLengthHilbert(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,filterOrder);
        length_temp_all_trials=[];
        for ii=1:length(length_temp_fg)
             length_temp_all_trials=[length_temp_all_trials,length_temp_fg{ii}];
        end
        length_gatherer_fg_hilbert{counter}=length_temp_all_trials;
        %%%%%% Wavelet transform %%%%%%%%%%%%%%%%%%%
        [length_temp_fg,~]= getBurstLengthWavelet(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
        length_temp_all_trials=[];
        for ii=1:length(length_temp_fg)
            length_temp_all_trials=[length_temp_all_trials,length_temp_fg{ii}];
        end
        length_gatherer_fg_wavelet{counter}=length_temp_all_trials;
        %%%%%% Feingold %%%%%%%%%%%%%%%%%%%
        [length_temp_fg]= getBurstLengthFeingold(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,filterOrder);
        length_temp_all_trials=[];
        for ii=1:length(length_temp_fg)
            length_temp_all_trials=[length_temp_all_trials,length_temp_fg{ii}];
        end
        length_gatherer_fg_feingold{counter}=length_temp_all_trials;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% OMP-GEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thresholdFactor=0.75*sqrt(diffPower);
        [length_temp_fg,~,~,~,~,~]= getBurstLength_all(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp,'OMP-GEAR');
        length_temp_all_trials=[];
        for ii=1:length(length_temp_fg)
            length_temp_all_trials=[length_temp_all_trials,length_temp_fg{ii}'];
        end
        length_gatherer_fg_OMP_gear{counter}=length_temp_all_trials;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        counter=counter+1;
    end
    length_all_elec_sg_hilbert=[];
    length_all_elec_fg_hilbert=[];
    median_gatherer_sg_hilbert=zeros(1,num_elec);
    median_gatherer_fg_hilbert=zeros(1,num_elec);
    length_all_elec_sg_wavelet=[];
    length_all_elec_fg_wavelet=[];
    median_gatherer_sg_wavelet=zeros(1,num_elec);
    median_gatherer_fg_wavelet=zeros(1,num_elec);
    length_all_elec_sg_feingold=[];
    length_all_elec_fg_feingold=[];
    median_gatherer_sg_feingold=zeros(1,num_elec);
    median_gatherer_fg_feingold=zeros(1,num_elec);
    length_all_elec_sg_OMP_gear=[];
    length_all_elec_fg_OMP_gear=[];
    median_gatherer_sg_OMP_gear=zeros(1,num_elec);
    median_gatherer_fg_OMP_gear=zeros(1,num_elec);
    for i=1:num_elec
        %%%%%%%%%%%%% Hilbert transform %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        reject_sg_indices=find(length_gatherer_sg_hilbert{i} > 0.8);
        reject_fg_indices=find(length_gatherer_fg_hilbert{i} > 0.8);
        length_gatherer_sg_hilbert{i}(reject_sg_indices)=[];
        length_gatherer_fg_hilbert{i}(reject_fg_indices)=[];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        median_gatherer_sg_hilbert(i)=median(length_gatherer_sg_hilbert{i});
        median_gatherer_fg_hilbert(i)=median(length_gatherer_fg_hilbert{i});
        if (length(length_gatherer_sg_hilbert{i})<30 || length(length_gatherer_fg_hilbert{i})<30)
            median_gatherer_sg_hilbert(i)=0;
            median_gatherer_fg_hilbert(i)=0;
            power_gatherer_sg_hilbert(i)=0;
            power_gatherer_fg_hilbert(i)=0;
        end
        %%%%%%%%%%%%% Wavelet transform %%%%%%%%%%%%%%%%%%%
        reject_sg_indices=find(length_gatherer_sg_wavelet{i} > 0.8);
        reject_fg_indices=find(length_gatherer_fg_wavelet{i} > 0.8);
        length_gatherer_sg_wavelet{i}(reject_sg_indices)=[];
        length_gatherer_fg_wavelet{i}(reject_fg_indices)=[];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        median_gatherer_sg_wavelet(i)=median(length_gatherer_sg_wavelet{i});
        median_gatherer_fg_wavelet(i)=median(length_gatherer_fg_wavelet{i});
        if length(length_gatherer_sg_wavelet{i})<20 || length(length_gatherer_fg_wavelet{i})<20
            median_gatherer_sg_wavelet(i)=0;
            median_gatherer_fg_wavelet(i)=0;
            power_gatherer_sg_wavelet(i)=0;
            power_gatherer_fg_wavelet(i)=0;
        end
        %%%%%%%%%%%%%% Feingold %%%%%%%%%%%%%%%%%%%%%
        reject_sg_indices=find(length_gatherer_sg_feingold{i} > 0.8);
        reject_fg_indices=find(length_gatherer_fg_feingold{i} > 0.8);
        length_gatherer_sg_feingold{i}(reject_sg_indices)=[];
        length_gatherer_fg_feingold{i}(reject_fg_indices)=[];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        median_gatherer_sg_feingold(i)=median(length_gatherer_sg_feingold{i});
        median_gatherer_fg_feingold(i)=median(length_gatherer_fg_feingold{i});
        if length(length_gatherer_sg_feingold{i})<20 || length(length_gatherer_fg_feingold{i})<20
            median_gatherer_sg_feingold(i)=0;
            median_gatherer_fg_feingold(i)=0;
            power_gatherer_sg_feingold(i)=0;
            power_gatherer_fg_feingold(i)=0;
        end
        %%%%% OMP-GEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        reject_sg_indices=find(length_gatherer_sg_OMP_gear{i} > 0.8);
        reject_fg_indices=find(length_gatherer_fg_OMP_gear{i} > 0.8);
        length_gatherer_sg_OMP_gear{i}(reject_sg_indices)=[];
        length_gatherer_fg_OMP_gear{i}(reject_fg_indices)=[];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        median_gatherer_sg_OMP_gear(i)=median(length_gatherer_sg_OMP_gear{i});
        median_gatherer_fg_OMP_gear(i)=median(length_gatherer_fg_OMP_gear{i});
        if (length(length_gatherer_sg_OMP_gear{i})<20) || (length(length_gatherer_fg_OMP_gear{i})<20)
            median_gatherer_sg_OMP_gear(i)=0;
            median_gatherer_fg_OMP_gear(i)=0;
            power_gatherer_sg_OMP_gear(i)=0;
            power_gatherer_fg_OMP_gear(i)=0;
        end
    end
    power_gatherer_sg_hilbert=power_gatherer_sg_hilbert(power_gatherer_sg_hilbert~=0);
    power_gatherer_fg_hilbert=power_gatherer_fg_hilbert(power_gatherer_fg_hilbert~=0);
    median_gatherer_sg_hilbert=median_gatherer_sg_hilbert(median_gatherer_sg_hilbert~=0);
    median_gatherer_fg_hilbert=median_gatherer_fg_hilbert(median_gatherer_fg_hilbert~=0);
    power_gatherer_sg_wavelet=power_gatherer_sg_wavelet(power_gatherer_sg_wavelet~=0);
    power_gatherer_fg_wavelet=power_gatherer_fg_wavelet(power_gatherer_fg_wavelet~=0);
    median_gatherer_sg_wavelet=median_gatherer_sg_wavelet(median_gatherer_sg_wavelet~=0);
    median_gatherer_fg_wavelet=median_gatherer_fg_wavelet(median_gatherer_fg_wavelet~=0);
    power_gatherer_sg_feingold=power_gatherer_sg_feingold(power_gatherer_sg_feingold~=0);
    power_gatherer_fg_feingold=power_gatherer_fg_feingold(power_gatherer_fg_feingold~=0);
    median_gatherer_sg_feingold=median_gatherer_sg_feingold(median_gatherer_sg_feingold~=0);
    median_gatherer_fg_feingold=median_gatherer_fg_feingold(median_gatherer_fg_feingold~=0);
    power_gatherer_sg_OMP_gear=power_gatherer_sg_OMP_gear(power_gatherer_sg_OMP_gear~=0);
    power_gatherer_fg_OMP_gear=power_gatherer_fg_OMP_gear(power_gatherer_fg_OMP_gear~=0);
    median_gatherer_sg_OMP_gear=median_gatherer_sg_OMP_gear(median_gatherer_sg_OMP_gear~=0);
    median_gatherer_fg_OMP_gear=median_gatherer_fg_OMP_gear(median_gatherer_fg_OMP_gear~=0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    power_gatherer_sg_hilbert=10*log10(power_gatherer_sg_hilbert);
    power_gatherer_fg_hilbert=10*log10(power_gatherer_fg_hilbert);
    power_gatherer_sg_wavelet=10*log10(power_gatherer_sg_wavelet);
    power_gatherer_fg_wavelet=10*log10(power_gatherer_fg_wavelet);
    power_gatherer_sg_feingold=10*log10(power_gatherer_sg_feingold);
    power_gatherer_fg_feingold=10*log10(power_gatherer_fg_feingold);
    power_gatherer_sg_OMP_gear=10*log10(power_gatherer_sg_OMP_gear);
    power_gatherer_fg_OMP_gear=10*log10(power_gatherer_fg_OMP_gear);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Wavelet transform %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,2))
    % [matched_sg_indices,matched_fg_indices]=power_matching_hist(power_gatherer_sg_wavelet,power_gatherer_fg_wavelet);
    matched_sg_lengths=median_gatherer_sg_wavelet(matched_sg_indices);
    matched_fg_lengths=median_gatherer_fg_wavelet(matched_fg_indices);
    violin_swarm_plot(matched_sg_lengths,matched_fg_lengths);
    if Monkey_num==1
        title('Wavelet','FontSize',24,'FontWeight','bold','FontName','Helvetica');
        save('length_wavelet_matched_power_M1.mat','matched_sg_lengths','matched_fg_lengths','matched_sg_indices','matched_fg_indices');
    end
    if Monkey_num==2
        save('length_wavelet_matched_power_M2.mat','matched_sg_lengths','matched_fg_lengths','matched_sg_indices','matched_fg_indices');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Feingold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
     subplot(plotHandles(Monkey_num,3))
    % [matched_sg_indices,matched_fg_indices]=power_matching_hist(power_gatherer_sg_feingold,power_gatherer_fg_feingold);
    matched_sg_lengths=median_gatherer_sg_feingold(matched_sg_indices);
    matched_fg_lengths=median_gatherer_fg_feingold(matched_fg_indices);
    violin_swarm_plot(matched_sg_lengths,matched_fg_lengths); 
%}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% OMP-GEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    matched_sg_lengths=median_gatherer_sg_OMP_gear(matched_sg_indices);
    matched_fg_lengths=median_gatherer_fg_OMP_gear(matched_fg_indices);
    violin_swarm_plot(matched_sg_lengths,matched_fg_lengths);
    if Monkey_num==1
        title('OMP-GEAR','FontSize',24,'FontWeight','bold','FontName','Helvetica');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
