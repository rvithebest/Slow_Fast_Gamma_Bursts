clc;clear
figure;
for Monkey_num=1:2
    clearvars -except Monkey_num
    % parent_file_path='C:\Users\VIgneshSRLab\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
    parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
    displayFlag=0;
    stimulusPeriodS=[0.25 0.75];
    baselinePeriodS=[-0.5 0];
    S_stim_accumulator=[];
    S_bl_accumulator=[];
    S_norm_accumulator=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Monkey_num==1
        % Monkey- alpaH
        load((fullfile(parent_file_path,'alpaH_info','parameterCombinations.mat')))
        load(fullfile(parent_file_path,'alpaH_info','badTrials.mat'));
        load(fullfile(parent_file_path,'alpaH_info','alpaHMicroelectrodeRFData.mat'));
        LFP_data_file=dir(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H'));
        LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'}));
        LFP_data_file = natsortfiles({LFP_data_file.name});
        current_electrode=highRMSElectrodes(1:77);
        select_pos=37;
        current_electrode=current_electrode(select_pos);
    elseif Monkey_num==2
        % Monkey- kesariH
        load((fullfile(parent_file_path,'kesariH_info','parameterCombinations.mat')))
        load(fullfile(parent_file_path,'kesariH_info','badTrials_kesari.mat'));
        load(fullfile(parent_file_path,'kesariH_info','kesariHMicroelectrodeRFData_Two.mat'));
        LFP_data_file=dir(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H'));
        LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'}));
        LFP_data_file = natsortfiles({LFP_data_file.name});
        current_electrode=highRMSElectrodes(1:31);
        % select_pos=9;
        % select_pos=2;
        select_pos=14;
        current_electrode=current_electrode(select_pos);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('timeVals.mat')
    num_elec=length(current_electrode);
    counter=1;
    for i=current_electrode
        if Monkey_num==1
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H',LFP_data_file{i}))
        end
        if Monkey_num==2
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H',LFP_data_file{i}))
        end
        trial_temp=[parameterCombinations{:,:,:,2,3},parameterCombinations{:,:,:,3,3},parameterCombinations{:,:,:,4,3}];
        % trial_temp=[parameterCombinations{:,:,:,2,8}];
        trial_temp=setdiff(trial_temp,badTrials);
        data_temp=analogDataDecimated(trial_temp,:);
        % FOR 100 Hz
        Fs=250;
        params.Fs=Fs;
        params.fpass=[0 100];
        params.tapers=[1 1];
        params.trialave=1;
        params.pad=-1;
        params.err=[2 0.05];
        s1_index=find(timeVals>stimulusPeriodS(1),1);
        s2_index=find(timeVals>stimulusPeriodS(2),1);
        LFP_signal_stim=data_temp(:,s1_index:s2_index);
        [S_stim,f_stim]=mtspectrumc(LFP_signal_stim',params);
        S_stim_accumulator=[S_stim_accumulator;S_stim'];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s1_index=find(timeVals>baselinePeriodS(1),1);
        s2_index=find(timeVals>baselinePeriodS(2),1);
        LFP_signal_bl=data_temp(:,s1_index:s2_index);
        [S_bl,f_bl]=mtspectrumc(LFP_signal_bl',params);
        S_bl_accumulator=[S_bl_accumulator;S_bl'];
        S_norm=S_stim./S_bl;
        
        S_norm_accumulator=[S_norm_accumulator;S_norm'];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        counter=counter+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_stim_avg=mean((S_stim_accumulator),1);
    hold on;
    S_bl_avg=mean((S_bl_accumulator),1);
    S_norm_avg=mean((S_norm_accumulator),1);
    % Take the average after substracting- correct power computation
    subplot(2,1,(1+(1*(Monkey_num-1))))
    plot(f_stim, 10*log10(S_norm_avg), '-k', 'LineWidth', 2);
    % plot(f_bl, 10*log10(S_stim_avg./S_bl_avg), '-k', 'LineWidth', 2);
    % plot(f_bl, log10(S_stim_avg), '-k', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    % set gca font size to 18
    set(gca,'FontSize',18);
    set(gca, 'FontName', 'Calibri');
    % Draw the frequency ranges of slow and fast gamma corresponding to M!
    % and M2
    if Monkey_num==1
        % Vertical line at slow gamma range
        % Slow gamma range for M1= [20 32]- Orange color
        orange_color=[1 0.6471 0];
        xline(20, 'Color', 'blue', 'LineWidth', 2);
        xline(32, 'Color', 'blue', 'LineWidth', 2);
        % Fast gamma range for M1= [36 65]- Blue color
        xline(36, 'Color', orange_color, 'LineWidth', 2);
        xline(65, 'Color', orange_color, 'LineWidth', 2);
    elseif Monkey_num==2
        orange_color=[1 0.6471 0];
        % Vertical line at slow gamma range
        % Slow gamma range for M2= [20 38]- Orange color
        xline(20, 'Color', 'blue', 'LineWidth', 2);
        xline(38, 'Color', 'blue', 'LineWidth', 2);
        % Fast gamma range for M2= [42 65]- Blue color
        xline(42, 'Color', orange_color, 'LineWidth', 2);
        xline(65, 'Color', orange_color, 'LineWidth', 2);
    end
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        annotation('textbox',...
    [0.05 0.68 0.12 0.12],...
    'String',{'Monkey 1'},...
    'Rotation',90,...
    'FontWeight','bold',...
    'FontSize',16,...
    'FontName','Helvetica',...
    'EdgeColor',[1 1 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     annotation('textbox',...
    [0.05 0.21 0.12 0.12],...
    'String',{'Monkey 2'},...
    'Rotation',90,...
    'FontWeight','bold',...
    'FontSize',16,...
    'FontName','Helvetica',...
    'EdgeColor',[1 1 1]);


























