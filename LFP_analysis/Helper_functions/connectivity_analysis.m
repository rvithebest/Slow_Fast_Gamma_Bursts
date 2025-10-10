clc;clear
f=figure;
f.WindowState="Maximized";
plotHandles=getPlotHandles(2,4,[0.08 0.08 0.9 0.9],0.07,0.06,0);
for Monkey_num=1:2
    clearvars -except Monkey_num f plotHandles
    parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
    displayFlag=0;
    stimulusPeriodS=[(0.25+(1/250)) 0.75];
    baselinePeriodS=[(-0.5+(1/250)) 0];
    baselinePeriodS=[(-0.5+(1/2000)) 0];
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
        LFP_data_file=dir(fullfile(parent_file_path,'LFP_data_alpa_H'));
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
        load(fullfile(parent_file_path,'kesariH_info','badTrials.mat'));
        load(fullfile(parent_file_path,'kesariH_info','kesariHMicroelectrodeRFData.mat'));
        LFP_data_file=dir(fullfile(parent_file_path,'LFP_data_kesari_H'));
        LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'}));
        LFP_data_file = natsortfiles({LFP_data_file.name});
        slow_gamma_freq=[20 38];
        fast_gamma_freq=[42 65];
        gabor_accumulator=gaborInfo_accumulator_kesariH;
        header_accumulator=header_accumulator_kesariH;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('timeVals.mat');
    num_elec=length(current_electrode);
    counter=1;
    SF_num=2;
    ORI_num=9; % all ORI
    trial_temp=parameterCombinations{:,:,:,SF_num,ORI_num};
    trial_temp = setdiff(trial_temp,badTrials);
    % Define a struct named data for FeildTrip compatibility
    data = struct();
    data.trial=cell(1,length(trial_temp));
    data.time=cell(1,length(trial_temp));
    data.label=cell(num_elec,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % data.fsample=250;
    data.fsample=2000;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=current_electrode
        if Monkey_num==1
            load(fullfile(parent_file_path,'LFP_data_alpa_H',LFP_data_file{i}))
        elseif Monkey_num==2
            load(fullfile(parent_file_path,'LFP_data_kesari_H',LFP_data_file{i}))
        end
        for j=1:length(trial_temp)
            data.trial{j}=[data.trial{j}; analogDataDecimated(trial_temp(j),:)];
            data.time{j}=timeVals;
        end
        data.label{counter}=['elec_'  num2str(i)];
        counter=counter+1;
    end
    %%%%%%%%%%%%%%%%%% Data preprocessing %%%%%%%%%%%%%%%%%%%%%%
    cfg=[];
    cfg.toilim=baselinePeriodS;
    data_bl=ft_redefinetrial(cfg,data);
    cfg=[];
    cfg.toilim=stimulusPeriodS;
    data_stim=ft_redefinetrial(cfg,data);
    %%%%%%%%%%%%%%%% PPC and WPLI %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfg=[];
    cfg.method='mtmfft';
    cfg.taper='dpss';
    cfg.tapsmofrq=4;
    anal_details=cfg;
    % Num tapers= (2*4*0.5)-1=3
    cfg.output='fourier';
    cfg.keeptrials='yes';
    cfg.foilim=[0 100];
    freq_bl=ft_freqanalysis(cfg,data_bl);
    freq_stim=ft_freqanalysis(cfg,data_stim);
    cfg=[];
    % WPLI or PLV
    cfg.method='wpli';
    % cfg.method='plv';
    anal_details_2=cfg;
    conn_stat_bl=ft_connectivityanalysis(cfg,freq_bl);
    conn_stat_stim=ft_connectivityanalysis(cfg,freq_stim);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Have to decide whether to take abs or not
    conn_bl=abs(conn_stat_bl.([cfg.method 'spctrm']));
    conn_stim=abs(conn_stat_stim.([cfg.method 'spctrm']));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Monkey_num==1
        save("M1_conn_analysis_WPLI.mat","conn_bl","conn_stim","freq_bl","freq_stim","anal_details","anal_details_2");
    elseif Monkey_num==2
        save("M2_conn_analysis_WPLI.mat","conn_bl","conn_stim","freq_bl","freq_stim","anal_details","anal_details_2");
    end
end