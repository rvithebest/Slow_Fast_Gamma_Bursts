clc;clear
f=figure;
f.WindowState="Maximized";
plotHandles=getPlotHandles(2,4,[0.08 0.08 0.9 0.9],0.07,0.06,0);
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
        gabor_accumulator=gaborInfo_accumulator_alpaH;
        header_accumulator=header_accumulator_alpaH;
        load('M1_power_gatherer.mat');
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
        gabor_accumulator=gaborInfo_accumulator_kesariH;
        header_accumulator=header_accumulator_kesariH;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('timeVals.mat')
    num_elec=length(current_electrode);
    counter=1;
    length_gatherer_sg=cell(1,num_elec);
    length_gatherer_fg=cell(1,num_elec);
    power_gatherer_sg=zeros(1,num_elec);
    power_gatherer_fg=zeros(1,num_elec);
    amp_gatherer_sg=cell(1,num_elec);
    amp_gatherer_fg=cell(1,num_elec);
    tf_accum_sg=[];
    tf_accum_fg=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        % power_gatherer_sg(counter)= get_Power(data_temp, timeVals,stimulusPeriodS, slow_gamma_freq);
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        tf_accum_sg=[tf_accum_sg,thresholdFactor];
        % thresholdFactor=sqrt(thresholdFraction);
        [length_temp_sg,~,~,~,~,a_temp_sg]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        length_temp_all_trials=[];
        a_temp_all_trials=[];
        for ii=1:length(length_temp_sg)
            length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}'];
            a_temp_all_trials=[a_temp_all_trials,a_temp_sg{ii}'];
        end
        length_gatherer_sg{counter}=length_temp_all_trials;
        amp_gatherer_sg{counter}=a_temp_all_trials;
        %%%%%% Fast gamma burst computation %%%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
        power_gatherer_fg(counter)=diffPower;
        % power_gatherer_fg(counter)= get_Power(data_temp, timeVals,stimulusPeriodS, fast_gamma_freq);
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        tf_accum_fg=[tf_accum_fg,thresholdFactor];
        % thresholdFactor=sqrt(thresholdFraction);
        [length_temp_fg,~,~,~,~,a_temp_fg]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        length_temp_all_trials=[];
        a_temp_all_trials=[]; 
        for ii=1:length(length_temp_fg)
            length_temp_all_trials=[length_temp_all_trials,length_temp_fg{ii}'];
            a_temp_all_trials=[a_temp_all_trials,a_temp_fg{ii}'];
        end
        length_gatherer_fg{counter}=length_temp_all_trials;
        amp_gatherer_fg{counter}=a_temp_all_trials;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        counter=counter+1;
    end
    length_all_elec_sg=[];
    length_all_elec_fg=[];
    median_gatherer_sg=zeros(1,num_elec);
    median_gatherer_fg=zeros(1,num_elec);
    median_gatherer_sg_amp=zeros(1,num_elec);
    median_gatherer_fg_amp=zeros(1,num_elec);
    amp_all_elec_sg=[];
    amp_all_elec_fg=[];
    for i=1:num_elec
        reject_idx=find(length_gatherer_sg{i}>0.8);
        amp_gatherer_sg{i}(reject_idx)=[];
        reject_idx=find(length_gatherer_fg{i}>0.8);
        amp_gatherer_fg{i}(reject_idx)=[];
        length_gatherer_sg{i}=length_gatherer_sg{i}(length_gatherer_sg{i}<0.8);
        length_gatherer_fg{i}=length_gatherer_fg{i}(length_gatherer_fg{i}<0.8);
        median_gatherer_sg(i)=median(length_gatherer_sg{i});
        median_gatherer_fg(i)=median(length_gatherer_fg{i});
        median_gatherer_sg_amp(i)=median(amp_gatherer_sg{i});
        median_gatherer_fg_amp(i)=median(amp_gatherer_fg{i});
        % if less than 20 bursts are detected for either slow or fast gamma, then exclude that electrode
        if ((length(length_gatherer_sg{i})<20) || (length(length_gatherer_fg{i})<20))
            median_gatherer_sg(i)=0;
            median_gatherer_fg(i)=0;
            power_gatherer_sg(i)=0;
            power_gatherer_fg(i)=0;
            median_gatherer_sg_amp(i)=0;
            median_gatherer_fg_amp(i)=0;
            continue;
        end
        length_all_elec_sg=[length_all_elec_sg,length_gatherer_sg{i}];
        length_all_elec_fg=[length_all_elec_fg,length_gatherer_fg{i}];
        amp_all_elec_sg=[amp_all_elec_sg,amp_gatherer_sg{i}];
        amp_all_elec_fg=[amp_all_elec_fg,amp_gatherer_fg{i}];
    end
    power_gatherer_sg=power_gatherer_sg(power_gatherer_sg~=0);
    power_gatherer_fg=power_gatherer_fg(power_gatherer_fg~=0);
    median_gatherer_sg=median_gatherer_sg(median_gatherer_sg~=0);
    median_gatherer_fg=median_gatherer_fg(median_gatherer_fg~=0);
    median_gatherer_sg_amp=median_gatherer_sg_amp(median_gatherer_sg_amp~=0);
    median_gatherer_fg_amp=median_gatherer_fg_amp(median_gatherer_fg_amp~=0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    power_gatherer_sg=10*log10(power_gatherer_sg);
    power_gatherer_fg=10*log10(power_gatherer_fg);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,1))
    histogram_all_burst(length_all_elec_sg,length_all_elec_fg,Monkey_num,10);
    if Monkey_num==1
        % save("MP_burst_duration_M1_fig_2.mat","length_all_elec_sg","length_all_elec_fg");
        % Perform a 2-sample Kolmogorov-Smirnov test to compare the distributions
        [h,pValue] = kstest2(length_all_elec_sg,length_all_elec_fg);
        disp(['Monkey 1: p-value = ', num2str(pValue),'for distribution comparison']);
    elseif Monkey_num==2
        % save("MP_burst_duration_M2_fig_2.mat","length_all_elec_sg","length_all_elec_fg");
        [h,pValue] = kstest2(length_all_elec_sg,length_all_elec_fg);
        disp(['Monkey 2: p-value = ', num2str(pValue),'for distribution comparison']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,2));
    % violin_swarm_plot(median_gatherer_sg,median_gatherer_fg);
    violin_swarm_plot_paired(median_gatherer_sg,median_gatherer_fg,0);
    if Monkey_num==1
       ylim([0 0.6])
    end 
    if Monkey_num==2
        ylim([0 0.4])
    end
    % title('B');
    if Monkey_num==1
        save('MP_results_M1_fig_2.mat','median_gatherer_sg','median_gatherer_fg','power_gatherer_sg','power_gatherer_fg');
    elseif Monkey_num==2
        save('MP_results_M2_fig_2.mat','median_gatherer_sg','median_gatherer_fg','power_gatherer_sg','power_gatherer_fg');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Power matching electrodes
    [matched_sg_indices,matched_fg_indices]=power_matching_hist(power_gatherer_sg,power_gatherer_fg);
    matched_sg_lengths = median_gatherer_sg(matched_sg_indices);
    matched_fg_lengths = median_gatherer_fg(matched_fg_indices);
    if Monkey_num==1
        save('matched_power_idx_M1_fig_2.mat','matched_sg_indices','matched_fg_indices','matched_sg_lengths','matched_fg_lengths');
    elseif Monkey_num==2
        save('matched_power_idx_M2_fig_2.mat','matched_sg_indices','matched_fg_indices','matched_sg_lengths','matched_fg_lengths');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,4));
    violin_swarm_plot(matched_sg_lengths,matched_fg_lengths)
    if Monkey_num==1
       ylim([0 0.6])
    end 
    if Monkey_num==2
        ylim([0 0.4])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,3))
    % scatter-plot of slow gamma, unfilled circles
    scatter_plot_power_burst(power_gatherer_sg,power_gatherer_fg,median_gatherer_sg ...
        ,median_gatherer_fg,matched_sg_indices,matched_fg_indices,Monkey_num)
    if Monkey_num==1
       ylim([0 0.6]);
    end 
    if Monkey_num==2
        ylim([0 0.4])
    end
    set_axis_ticks_fontsize(plotHandles,22,18,Monkey_num);
end
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