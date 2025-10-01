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
    load('timeVals.mat')
    num_elec=length(current_electrode);
    counter=1;
    length_gatherer_sg=cell(8,num_elec);
    length_gatherer_fg=cell(8,num_elec);
    power_gatherer_sg_ORI=zeros(8,num_elec);
    power_gatherer_fg_ORI=zeros(8,num_elec);
    for i=current_electrode
        if Monkey_num==1
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H',LFP_data_file{i}))
        end
        if Monkey_num==2
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H',LFP_data_file{i}))
        end
        power_stim_sg=zeros(1,8);power_bl_sg=zeros(1,8);
        power_stim_fg=zeros(1,8);power_bl_fg=zeros(1,8);
        SF_num=2;
        for ORI=1:8
            trial_temp=parameterCombinations{:,:,:,SF_num,ORI};
            trial_temp=setdiff(trial_temp,badTrials);
            data_temp=analogDataDecimated(trial_temp,:);
            gabor_temp=gabor_accumulator{SF_num,ORI,counter};
            header_temp=header_accumulator{SF_num,ORI,counter};
            %%%%%% Power computation %%%%%%
            power_stim_sg(1,ORI)=log10(get_Power(data_temp,timeVals,stimulusPeriodS,slow_gamma_freq));
            power_bl_sg(1,ORI)=log10(get_Power(data_temp,timeVals,baselinePeriodS,slow_gamma_freq));
            power_stim_fg(1,ORI)=log10(get_Power(data_temp,timeVals,stimulusPeriodS,fast_gamma_freq));
            power_bl_fg(1,ORI)=log10(get_Power(data_temp,timeVals,baselinePeriodS,fast_gamma_freq));
            %%%%%%%% Slow gamma burst computation %%%%%%
            diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,slow_gamma_freq);
            thresholdFactor=thresholdFraction*sqrt(diffPower);
            [length_temp_sg,~,~,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
            length_temp_all_trials=[];
            for ii=1:length(length_temp_sg)
                length_temp_all_trials=[length_temp_all_trials,length_temp_sg{ii}'];
            end
            length_gatherer_sg{ORI,counter}=length_temp_all_trials;
            %%%%%%%% Fast gamma burst computation %%%%%%
            diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
            thresholdFactor=thresholdFraction*sqrt(diffPower);
            [length_temp_fg,~,~,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
            length_temp_all_trials=[];
            for ii=1:length(length_temp_fg)
                length_temp_all_trials=[length_temp_all_trials,length_temp_fg{ii}'];
            end
            length_gatherer_fg{ORI,counter}=length_temp_all_trials;
        end
        bl_avg_sg=mean(power_bl_sg,2);
        bl_avg_fg=mean(power_bl_fg,2);
        for ORI=1:8
            power_gatherer_sg_ORI(ORI,counter)=10*(power_stim_sg(1,ORI)-bl_avg_sg);
            power_gatherer_fg_ORI(ORI,counter)=10*(power_stim_fg(1,ORI)-bl_avg_fg);
        end
        counter=counter+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Monkey_num==1
         save("M1_power_gatherer.mat","power_gatherer_sg_ORI","power_gatherer_fg_ORI");
    elseif Monkey_num==2
         save("M2_power_gatherer.mat","power_gatherer_sg_ORI","power_gatherer_fg_ORI");
    end
    length_all_elec_sg=[];
    length_all_elec_fg=[];
    median_gatherer_sg=zeros(8,num_elec);
    median_gatherer_fg=zeros(8,num_elec);
    for ii=1:num_elec
        for ORI=1:8
            length_gatherer_sg{ORI,ii}=length_gatherer_sg{ORI,ii}(length_gatherer_sg{ORI,ii}<0.8);
            length_gatherer_fg{ORI,ii}=length_gatherer_fg{ORI,ii}(length_gatherer_fg{ORI,ii}<0.8);
            median_gatherer_sg(ORI,ii)=median(length_gatherer_sg{ORI,ii});
            median_gatherer_fg(ORI,ii)=median(length_gatherer_fg{ORI,ii});
            % if less than 20 bursts are detected for either slow or fast gamma, then exclude that electrode
            if (length(length_gatherer_sg{ORI,ii})<20) || (length(length_gatherer_fg{ORI,ii})<20)
                median_gatherer_sg(ORI,ii)=0;
                median_gatherer_fg(ORI,ii)=0;
                power_gatherer_sg_ORI(ORI,ii)=0;
                power_gatherer_fg_ORI(ORI,ii)=0;
                continue;
            end
            length_all_elec_sg=[length_all_elec_sg,length_gatherer_sg{ORI,ii}]; 
            length_all_elec_fg=[length_all_elec_fg,length_gatherer_fg{ORI,ii}];
        end
    end
    avg_power_sg=zeros(8,1); avg_power_fg=zeros(8,1);
    SEM_power_sg=zeros(8,1); SEM_power_fg=zeros(8,1);
    median_median_sg=zeros(8,1); median_median_fg=zeros(8,1);
    SEM_median_sg=zeros(8,1); SEM_median_fg=zeros(8,1);
    p_val=zeros(8,1);
    h_val=zeros(8,1);
    stats_val=cell(8,1);
    N_val=zeros(8,1);
    for ORI=1:8
        % take all electrodes with non-zero median burst length
        temp_length_sg=median_gatherer_sg(ORI,:);
        temp_length_sg=temp_length_sg(temp_length_sg~=0);
        temp_length_fg=median_gatherer_fg(ORI,:);
        temp_length_fg=temp_length_fg(temp_length_fg~=0);
        % 2-sample unpaired t-test
        % [h,p,ci,stats]=ttest2(temp_length_sg,temp_length_fg,'VarType','unequal');
        % 2-sample paired t-test
        % [h,p,ci,stats]=ttest(temp_length_sg,temp_length_fg,'Alpha',0.05,'Tail','both');
        N_temp=length(temp_length_fg);
        N_val(ORI)=N_temp;
        [p,h,stats] = signrank(temp_length_sg, temp_length_fg, 'alpha',0.05,'tail','right');
        p_val(ORI)=p;
        h_val(ORI)=h;
        stats_val{ORI}=stats;
        disp(['p-value for ORI ',num2str(ORI),' is ',num2str(p),' for monkey ',num2str(Monkey_num)]);
        temp_power_sg=power_gatherer_sg_ORI(ORI,:);
        temp_power_sg=temp_power_sg(temp_power_sg~=0);
        temp_power_fg=power_gatherer_fg_ORI(ORI,:);
        temp_power_fg=temp_power_fg(temp_power_fg~=0);
        avg_power_sg(ORI)=mean(temp_power_sg);
        avg_power_fg(ORI)=mean(temp_power_fg);
        SEM_power_sg(ORI)=std(temp_power_sg)/sqrt(length(temp_power_sg));
        SEM_power_fg(ORI)=std(temp_power_fg)/sqrt(length(temp_power_fg));
        median_median_sg(ORI)=median(temp_length_sg);
        median_median_fg(ORI)=median(temp_length_fg);
        SEM_median_sg(ORI)=getSEMedian(temp_length_sg,1000);
        SEM_median_fg(ORI)=getSEMedian(temp_length_fg,1000);
    end
    if Monkey_num==1
       save("stats_fig_3_M1.mat","p_val","h_val","stats_val","N_val");
    elseif Monkey_num==2
       save("stats_fig_3_M2.mat","p_val","h_val","stats_val","N_val");
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,1));
    histogram_all_burst(length_all_elec_sg,length_all_elec_fg,Monkey_num,10);
    % 2-sample Kolmogorov-Smirnov test to compare the distributions of burst lengths
    [h_ks,p_ks,ks_stat]=kstest2(length_all_elec_sg,length_all_elec_fg);
    disp(['Kolmogorov-Smirnov test p-value for monkey ',num2str(Monkey_num),' is ',num2str(p_ks),'for distribution of burst lengths']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ORI_axis=[0,22.5,45,67.5,90,112.5,135,157.5];
    ORI_axis=categorical(ORI_axis);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,3));
    errorbar(ORI_axis,avg_power_sg,SEM_power_sg,'-o','MarkerSize',6,'LineWidth',1.5,'Capsize',10,'Color','b');
    color_orange=[1 0.5 0];
    hold on;
    errorbar(ORI_axis,avg_power_fg,SEM_power_fg,'-o','MarkerSize',6,'LineWidth',1.5,'Capsize',10,'Color',color_orange);
    ylabel('Power (dB)');
    xlabel('Orientation (\circ)');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(Monkey_num,2));
    % median of burst lengths (median of medians)- for each ORI
    errorbar(ORI_axis,median_median_sg,SEM_median_sg,'-o','MarkerSize',6,'LineWidth',1.5,'Capsize',10,'Color','b');
    hold on;
    errorbar(ORI_axis,median_median_fg,SEM_median_fg,'-o','MarkerSize',6,'LineWidth',1.5,'Capsize',10,'Color',color_orange);
    if Monkey_num==1
        ylim([0 0.6]);
    elseif Monkey_num==2
        ylim([0 0.4]);
    end
    signific_val=p_val<0.05;
    %%% For significant p-values, the markerface color of that corresponding ORI is changed to red
    for i=1:8
        if signific_val(i)==1
            plot(ORI_axis(i),median_median_sg(i),'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','b','LineWidth',1.5,'Color','b');
            plot(ORI_axis(i),median_median_fg(i),'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor',color_orange,'LineWidth',1.5,'Color',color_orange);
        end
    end
    ylabel('Burst length (s)');
    xlabel('Orientation (\circ)');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set_axis_ticks_fontsize(plotHandles,22,16,Monkey_num);
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