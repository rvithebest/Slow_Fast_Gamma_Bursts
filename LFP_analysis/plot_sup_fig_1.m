clc;clear
f=figure;
f.WindowState="Maximized";
plotHandles_1=getPlotHandles(2,1,[0.08 0.08 0.15 0.86],0.07,0.06,0);
plotHandles_2=getPlotHandles(2,4,[0.3 0.54 0.65 0.4],0.02,0.04,0);
plotHandles_3=getPlotHandles(2,4,[0.3 0.08 0.65 0.4],0.02,0.04,0);
for Monkey_num=1:2
    clearvars -except Monkey_num plotHandles_1 plotHandles_2 plotHandles_3
    % parent_file_path='C:\Users\VigneshSRay\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
    parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
    displayFlag=0;
    stimulusPeriodS=[0.25 0.75];
    baselinePeriodS=[-0.5 0];
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
        select_pos=[37];
        % select_pos=54;
        current_electrode_selec=current_electrode(select_pos);
        plotHandle_temp=plotHandles_2;
        slow_gamma_freq=[20 32]; % M1
        fast_gamma_freq=[36 65]; % M1
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
        % select_pos=2;
        % select_pos=3,7;- new_Select
        select_pos=31;
        current_electrode_selec=current_electrode(select_pos);
        plotHandle_temp=plotHandles_3;
        slow_gamma_freq=[20 38]; % M2
        fast_gamma_freq=[42 65]; % M2
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('timeVals.mat')
    num_elec=length(current_electrode);
    counter=1;
    % FOR 100 Hz
    Fs=250;
    params.Fs=Fs;
    params.fpass=[0 100];
    params.tapers=[1 1];
    params.trialave=1;
    params.pad=-1;
    params.err=[2 0.05];
    for i=current_electrode_selec
        if Monkey_num==1
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H',LFP_data_file{i}))
        end
        if Monkey_num==2
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H',LFP_data_file{i}))
        end
        % trial_temp=[parameterCombinations{:,:,:,2,3},parameterCombinations{:,:,:,2,6},parameterCombinations{:,:,:,2,8}];
        trial_temp=[parameterCombinations{:,:,:,2,9}]; % Across all ORI
        trial_temp=setdiff(trial_temp,badTrials);
        data_temp=analogDataDecimated(trial_temp,:);
        s1_index=find(timeVals>stimulusPeriodS(1),1);
        s2_index=find(timeVals>stimulusPeriodS(2),1);
        LFP_signal_stim=data_temp(:,s1_index:s2_index);
        [S_stim,f_stim]=mtspectrumc(LFP_signal_stim',params);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s1_index=find(timeVals>baselinePeriodS(1),1);
        s2_index=find(timeVals>baselinePeriodS(2),1);
        LFP_signal_bl=data_temp(:,s1_index:s2_index);
        [S_bl,f_bl]=mtspectrumc(LFP_signal_bl',params);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        counter=counter+1;
    end
    movingWin=[0.25 0.025];
    S_accumulator=cell(num_elec,8);
    t_accumulator=cell(num_elec,8);
    f_accumulator=cell(num_elec,8);
    index=1;
    for i=current_electrode
        if Monkey_num==1
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H',LFP_data_file{i}))
        end
        if Monkey_num==2
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H',LFP_data_file{i}))
        end
        for ORI=1:8
            trial_temp=[parameterCombinations{:,:,:,2,ORI}]; 
            trial_temp=setdiff(trial_temp,badTrials);
            data_temp=analogDataDecimated(trial_temp,:);
            [S_temp,t_temp,f_temp]= mtspecgramc(data_temp',movingWin,params);
            t_temp = t_temp+timeVals(1)-(1/Fs); % Center the times with respect to the stimulus onset time
            [S_dB]=normalize_TF(S_temp,t_temp,baselinePeriodS);
            S_accumulator{index,ORI}=S_dB;
            t_accumulator{index,ORI}=t_temp;
            f_accumulator{index,ORI}=f_temp;
        end 
        index=index+1;
    end 
    S_trial_avg=cell(1,8);
    % ORI_axis=[ 0^{o},22.5^{o},45^{o},67.5^{o},90^{o},112.5^{o},135^{o},157.5^{o}];
    ORI_axis=[0,22.5,45,67.5,90,112.5,135,157.5];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ORI=1:8
        S_temp=cat(3,S_accumulator{:,ORI});
        t_temp=cell2mat(t_accumulator(1,ORI));
        f_temp=cell2mat(f_accumulator(1,ORI));
        S_trial_avg{1,ORI}=mean(S_temp,3);
        %%% TF plots %%%
        if ORI<5
            subplot(plotHandle_temp(1,ORI));
        else
            subplot(plotHandle_temp(2,ORI-4));
            if Monkey_num==2
                xlabel('Time (s)');
            end
        end
        pcolor(t_temp,f_temp,S_trial_avg{1,ORI});
        shading interp;
        colormap('jet');
        x=min(S_trial_avg{1,ORI}(:));
        y=max(S_trial_avg{1,ORI}(:));
        % caxis([x y-0.5]);
        disp(["M1 ORI "+num2str(ORI)+" min value: "+num2str(x)+" max value: "+num2str(y)]);
        if Monkey_num==1
           caxis([-3 8]);
        else 
           caxis([-4 10]);
        end 
        xlim([-0.1 1])
        ylim([0 100])
        if ORI<5 
            % Remove Xticks
            set(gca, 'XTick', []);
        else 
            if Monkey_num==1
                % Remove Xtick labels
                set(gca, 'XTick', []);
            end
            if Monkey_num==2
                xlabel('Time (s)');
                % Xticks- [0 0.8]
                set(gca, 'XTick', [0 0.8]);
            end
        end
        % draw black line at 0.25 s and 0.75 s (vertical lines)
        xline(0.25, 'Color', 'k', 'LineWidth', 2);
        xline(0.75, 'Color', 'k', 'LineWidth', 2);
        % draw horizontal line for slow_gamma_freq and fast_gamma_freq-black
        yline(slow_gamma_freq(1), 'Color', 'black', 'LineWidth', 2);
        yline(slow_gamma_freq(2), 'Color', 'black', 'LineWidth', 2);
        yline(fast_gamma_freq(1), 'Color', 'black', 'LineWidth', 2);
        yline(fast_gamma_freq(2), 'Color', 'black', 'LineWidth', 2);
        % Remove Ytick labels
        set(gca, 'YTick', []);
        if ORI==1
            ylabel('Frequency (Hz)');
            % Yticks- 0,50,100
            set(gca, 'YTick', [0 50 100]);
        end
        if ORI==5
            % Yticks- 0,50,100
            set(gca, 'YTick', [0 50 100]);
        end
        title([num2str(ORI_axis(ORI)),'^{o}'],"Interpreter",'tex','FontSize',16,'FontWeight','bold');
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Take the average after substracting- correct power computation
    subplot(plotHandles_1(Monkey_num,1));
    % plot(f_stim, 10*log10(S_norm_avg), '-k', 'LineWidth', 2);
    plot(f_bl, 10*log10(S_stim./S_bl), '-k', 'LineWidth', 2);
    % plot(f_bl, log10(S_stim_avg), '-k', 'LineWidth', 2);
    ylabel('Power (dB)');
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
        xlabel('Frequency (Hz)');
    end
    %%%%%%%%%%%% Time frequency plots (average) %%%%%%%%%%%%%%
    set_axis_ticks_fontsize(plotHandles_1,22,18,Monkey_num);
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set_axis_ticks_fontsize(plotHandles_2,22,18,1);
    set_axis_ticks_fontsize(plotHandles_2,22,18,2);
    set_axis_ticks_fontsize(plotHandles_3,22,18,1);
    set_axis_ticks_fontsize(plotHandles_3,22,18,2);
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
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     labels = {'A','B','C','D'};
     x_positions = [0.04, 0.25];
     y_positions = [0.9, 0.4];  % top and bottom rows
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


























