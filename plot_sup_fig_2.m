clc;clear
figure;
for Monkey_num=1:2
    clearvars -except Monkey_num
    parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
    displayFlag=0;
    stimulusPeriodS=[0.25 0.75];
    baselinePeriodS=[-0.5 0];
    thresholdFraction=0.5;
    num_iterations=120;
    dict_size=2500000;
    bin_num=10;
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
    load(LFP_data_file{91});
    timeVals=timeVals_decimated;
    num_elec=length(current_electrode);
    counter=1;
    SF_num=5;
    length_gatherer_sg=cell(SF_num,num_elec);
    length_gatherer_fg=cell(SF_num,num_elec);
    power_gatherer_sg_SF=zeros(SF_num,num_elec);
    power_gatherer_fg_SF=zeros(SF_num,num_elec);
    for i=current_electrode
        if Monkey_num==1
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H',LFP_data_file{i}))
        end
        if Monkey_num==2
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H',LFP_data_file{i}))
        end
        power_stim_sg=zeros(1,SF_num);power_bl_sg=zeros(1,SF_num);
        power_stim_fg=zeros(1,SF_num);power_bl_fg=zeros(1,SF_num);
        % At 157.5 degree
        ORI_num=8;
        for SF=1:5
            trial_temp=parameterCombinations{:,:,:,SF,ORI_num};
            trial_temp=setdiff(trial_temp,badTrials);
            data_temp=analogDataDecimated(trial_temp,:);
            gabor_temp=gabor_accumulator{SF,ORI_num,counter};
            header_temp=header_accumulator{SF,ORI_num,counter};
            %%%%%% Power computation %%%%%%
            power_stim_sg(1,SF)=log10(get_Power(data_temp,timeVals,stimulusPeriodS,slow_gamma_freq));
            power_bl_sg(1,SF)=log10(get_Power(data_temp,timeVals,baselinePeriodS,slow_gamma_freq));
            power_stim_fg(1,SF)=log10(get_Power(data_temp,timeVals,stimulusPeriodS,fast_gamma_freq));
            power_bl_fg(1,SF)=log10(get_Power(data_temp,timeVals,baselinePeriodS,fast_gamma_freq));
            %%%%%%%% Slow gamma burst computation %%%%%%
            diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,slow_gamma_freq);
            thresholdFactor=sqrt(thresholdFraction*diffPower);
            [length_temp_sg,~,~,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
            length_temp_all_trials=[];
            for ii=1:length(length_temp_sg)
                length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}'];
            end
            length_gatherer_sg{SF,counter}=length_temp_all_trials;
            %%%%%%%% Fast gamma burst computation %%%%%%%%%
            diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
            thresholdFactor=sqrt(thresholdFraction*diffPower);
            [length_temp_fg,~,~,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
            length_temp_all_trials=[];
            for ii=1:length(length_temp_fg)
                length_temp_all_trials=[length_temp_all_trials,length_temp_fg{ii}'];
            end
            length_gatherer_fg{SF,counter}=length_temp_all_trials;
        end
        bl_avg_sg=mean(power_bl_sg,2);
        bl_avg_fg=mean(power_bl_fg,2);
        for SF=1:5
            power_gatherer_sg_SF(SF,counter)=10*(power_stim_sg(1,SF)-bl_avg_sg);
            power_gatherer_fg_SF(SF,counter)=10*(power_stim_fg(1,SF)-bl_avg_fg);
        end
        counter=counter+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    length_all_elec_sg=[];
    length_all_elec_fg=[];
    median_gatherer_sg=zeros(5,num_elec);
    median_gatherer_fg=zeros(5,num_elec);
    for ii=1:num_elec
        for SF=1:5
            length_gatherer_sg{SF,ii}=length_gatherer_sg{SF,ii}(length_gatherer_sg{SF,ii}<0.8);
            length_gatherer_fg{SF,ii}=length_gatherer_fg{SF,ii}(length_gatherer_fg{SF,ii}<0.8);
            median_gatherer_sg(SF,ii)=median(length_gatherer_sg{SF,ii});
            median_gatherer_fg(SF,ii)=median(length_gatherer_fg{SF,ii});
            % if less than 20 bursts are detected for either slow or fast gamma, then exclude that electrode
            if length(length_gatherer_sg{SF,ii})<20 || length(length_gatherer_fg{SF,ii})<20
                median_gatherer_sg(SF,ii)=0;
                median_gatherer_fg(SF,ii)=0;
                power_gatherer_sg_SF(SF,ii)=0;
                power_gatherer_fg_SF(SF,ii)=0;
                continue;
            end
            length_all_elec_sg=[length_all_elec_sg,length_gatherer_sg{SF,ii}]; 
            length_all_elec_fg=[length_all_elec_fg,length_gatherer_fg{SF,ii}];
        end
    end
    avg_power_sg=zeros(5,1); avg_power_fg=zeros(5,1);
    SEM_power_sg=zeros(5,1); SEM_power_fg=zeros(5,1);
    median_median_sg=zeros(5,1); median_median_fg=zeros(5,1);
    SEM_median_sg=zeros(5,1); SEM_median_fg=zeros(5,1);
    p_val=zeros(5,1);
    median_matched_length_sg=zeros(5,1);
    median_matched_length_fg=zeros(5,1);
    SEM_matched_length_sg=zeros(5,1);
    SEM_matched_length_fg=zeros(5,1);
    p_val_2=zeros(5,1);
    for SF=1:5
        % take all electrodes with non-zero median burst length
        temp_length_sg=median_gatherer_sg(SF,:);
        temp_length_sg=temp_length_sg(temp_length_sg~=0);
        temp_length_fg=median_gatherer_fg(SF,:);
        temp_length_fg=temp_length_fg(temp_length_fg~=0);
        % 2-smaple unpaired t-test
        [h,p,ci,stats]=ttest2(temp_length_sg,temp_length_fg,"Vartype","unequal");
        p_val(SF)=p;
        disp(["p-value for SF " + num2str(SF) + ": " + num2str(p),"for- monkey " + num2str(Monkey_num)]);
        temp_power_sg=power_gatherer_sg_SF(SF,:);
        temp_power_sg=temp_power_sg(temp_power_sg~=0);
        temp_power_fg= power_gatherer_fg_SF(SF,:);
        temp_power_fg=temp_power_fg(temp_power_fg~=0);
        avg_power_sg(SF)=mean(temp_power_sg);
        avg_power_fg(SF)=mean(temp_power_fg);
        SEM_power_sg(SF)=std(temp_power_sg)/sqrt(length(temp_power_sg));
        SEM_power_fg(SF)=std(temp_power_fg)/sqrt(length(temp_power_fg));
        median_median_sg(SF)=median(temp_length_sg);
        median_median_fg(SF)=median(temp_length_fg);
        SEM_median_sg(SF)=getSEMedian(temp_length_sg);
        SEM_median_fg(SF)=getSEMedian(temp_length_fg);
        % Power matching does not lead to sufficient number of remaining
        % electrodes
        % Try power matching once more
        % [matched_sg_indices,matched_fg_indices]=power_matching_hist(temp_power_sg,temp_power_fg);
        % matched_length_sg= temp_length_sg(matched_sg_indices);
        % matched_length_fg= temp_length_fg(matched_fg_indices);
        % % 2-smaple unpaired t-test
        % [h,p,ci,stats]=ttest2(matched_length_sg,matched_length_fg,"Vartype","unequal");
        % p_val_2(SF)=p;
        % disp(["p-value for power matched SF " + num2str(SF) + ": " + num2str(p),"for- monkey " + num2str(Monkey_num)]);
        % median_matched_length_sg(SF)=median(matched_length_sg);
        % median_matched_length_fg(SF)=median(matched_length_fg);
        % SEM_matched_length_sg(SF)=getSEMedian(matched_length_sg);
        % SEM_matched_length_fg(SF)=getSEMedian(matched_length_fg);

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SF_axis=[0.5,1,2,4,8];
    SF_axis=categorical(SF_axis);
    subplot(2,3,(1+(3*(Monkey_num-1))))
    histogram_all_burst(length_all_elec_sg,length_all_elec_fg,Monkey_num,bin_num);
    % set gca font size to 18
    % errorbar(SF_axis,median_matched_length_sg,SEM_matched_length_sg,'-o','MarkerSize',6,'LineWidth',1.5,'Capsize',10,'Color','b');
    % hold on;
    % errorbar(SF_axis,median_matched_length_fg,SEM_matched_length_fg,'-o','MarkerSize',6,'LineWidth',1.5,'Capsize',10,'Color',[1 0.5 0]);
    set(gca,'FontSize',18);
    % set(gca, 'FontName', 'Calibri');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,3,(3+(3*(Monkey_num-1))))
    errorbar(SF_axis,avg_power_sg,SEM_power_sg,'-o','MarkerSize',6,'LineWidth',1.5,'Capsize',10,'Color','b');
    color_orange=[1 0.5 0];
    hold on;
    errorbar(SF_axis,avg_power_fg,SEM_power_fg,'-o','MarkerSize',6,'LineWidth',1.5,'Capsize',10,'Color',color_orange);
    ylabel('Power (dB)');
    xlabel('SF (cpd)');
    % set gca font size to 18
    set(gca,'FontSize',18);
    set(gca, 'FontName', 'Calibri');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,3,(2+(3*(Monkey_num-1))))
    % median of burst lengths (median of medians)- for each SF
    errorbar(SF_axis,median_median_sg,SEM_median_sg,'-o','MarkerSize',6,'LineWidth',1.5,'Capsize',10,'Color','b');
    hold on;
    errorbar(SF_axis,median_median_fg,SEM_median_fg,'-o','MarkerSize',6,'LineWidth',1.5,'Capsize',10,'Color',color_orange);
    ylabel('Burst length (s)');
    title('C');
    % set gca font size to 18
    set(gca,'FontSize',18);
    set(gca, 'FontName', 'Calibri');
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