function plot_bursts_PSD(plotHandles,Monkey_num)
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
    PSD_gatherer_sg_burst=zeros(size(freq_gatherer_sg_burst));
    PSD_gatherer_fg_burst=zeros(size(freq_gatherer_fg_burst));
    random_PSD_gatherer_sg_burst=zeros(size(freq_gatherer_sg_burst));
    random_PSD_gatherer_fg_burst=zeros(size(freq_gatherer_fg_burst));
    freq_gatherer_sg_burst=[];
    freq_gatherer_fg_burst=[];
    E_TF_MP_gatherer=cell(1,num_elec);
    E_avg_MP_gatherer=cell(1,num_elec);
    trial_idx_sg=1;
    trial_idx_fg=1;
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
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        [length_temp_sg,freq_temp_sg,time_center_temp_sg,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        [meanE,freqVals]=getEnergyMP3p1(gabor_temp,header_temp,timeVals);
        E_avg_MP_gatherer{counter}=meanE;
        meanE=E_avg_MP_gatherer{counter};
        threshold_pow=10^(-16);
        meanE(meanE<threshold_pow)=threshold_pow;
        meanE=log10(meanE);
        bl_idx_1=find(timeVals>=baselinePeriodS(1),1);
        bl_idx_2=find(timeVals<baselinePeriodS(2),1,'last');
        bl_avg=mean(meanE(:,bl_idx_1:bl_idx_2),2);
        E_TF_MP_gatherer{counter}=cell(1,length(length_temp_sg));
        for ii=1:length(length_temp_sg)
            E=mp2tf(squeeze(gabor_temp(ii,:,:)),header_temp(ii,:));
            % E= 256 (freq vals) x 512 (time vals)
            E_TF_MP_gatherer{counter}{ii}=E;
            E(E<threshold_pow)=threshold_pow;
            if isempty(length_temp_sg{ii}')
                continue;
            end
            reject_idx=find((length_temp_sg{ii}')>0.8);
            length_temp_sg{ii}(reject_idx)=[];
            time_center_temp_sg{ii}(reject_idx)=[];
            freq_temp_sg{ii}(reject_idx)=[];
            num_bursts=length(time_center_temp_sg{ii});
            for jj=1:num_bursts
                [~,time_idx]=min(abs(timeVals-time_center_temp_sg{ii}(jj)));
                PSD_gatherer_sg_burst(:,trial_idx_sg)=10*(log10(E(:,time_idx))-bl_avg);
                % Generate a random time point within the stimulus period
                random_time_center=stimulusPeriodS(1) + (stimulusPeriodS(2) - stimulusPeriodS(1)) * rand;
                [~,random_time_idx]=min(abs(timeVals-random_time_center));
                random_PSD_gatherer_sg_burst(:,trial_idx_sg)=10*(log10(E(:,random_time_idx))-bl_avg);
                trial_idx_sg=trial_idx_sg+1;
                freq_gatherer_sg_burst=[freq_gatherer_sg_burst,freqVals'];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Fast gamma burst computation %%%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        [length_temp_fg,freq_temp_fg,time_center_temp_fg,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        [~,freqVals]=getEnergyMP3p1(gabor_temp,header_temp,timeVals);
        for ii=1:length(length_temp_fg)
            if isempty(length_temp_fg{ii})
                continue;
            end
            E=mp2tf(squeeze(gabor_temp(ii,:,:)),header_temp(ii,:));
            % E= 256 (freq vals) x 512 (time vals)
            E=E_TF_MP_gatherer{counter}{ii};
            E(E<threshold_pow)=threshold_pow;
            reject_idx=find((length_temp_fg{ii}')>0.8);
            length_temp_fg{ii}(reject_idx)=[];
            time_center_temp_fg{ii}(reject_idx)=[];
            freq_temp_fg{ii}(reject_idx)=[];
            num_bursts=length(time_center_temp_fg{ii});
            for jj=1:num_bursts
                [~,time_idx]=min(abs(timeVals-time_center_temp_fg{ii}(jj)));
                freq_gatherer_fg_burst=[freq_gatherer_fg_burst,freqVals'];
                PSD_gatherer_fg_burst(:,trial_idx_fg)=10*(log10(E(:,time_idx))-bl_avg);
                % Generate a random time point within the stimulus period
                random_time_center=stimulusPeriodS(1) + (stimulusPeriodS(2) - stimulusPeriodS(1)) * rand;
                [~,random_time_idx]=min(abs(timeVals-random_time_center));
                random_PSD_gatherer_fg_burst(:,trial_idx_fg)=10*(log10(E(:,random_time_idx))-bl_avg);
                trial_idx_fg=trial_idx_fg+1;
            end
        end
        counter=counter+1;
        disp(counter);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Already computed and saved the TF values for all bursts in both monkeys. 
    if Monkey_num==1
         save('MP_TF_vals_M1.mat','E_TF_MP_gatherer','E_avg_MP_gatherer','freqVals','freq_gatherer_fg_burst','freq_gatherer_sg_burst','-v7.3');
    elseif Monkey_num==2
         save('MP_TF_vals_M2.mat','E_TF_MP_gatherer','E_avg_MP_gatherer','freqVals','freq_gatherer_fg_burst','freq_gatherer_sg_burst','-v7.3');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avg_window=9; % 4 Hz (freq resolution-0.5 Hz)
    subplot(plotHandles(Monkey_num,2));
    avg_PSD_gatherer_sg_burst=mean(PSD_gatherer_sg_burst,2);
    avg_PSD_gatherer_sg_burst=movmean(avg_PSD_gatherer_sg_burst,avg_window);
    plot(freqVals,avg_PSD_gatherer_sg_burst,'LineWidth',2,'Color','b');
    hold on;
    color_orange=[1,0.5,0];
    avg_PSD_gatherer_fg_burst=mean(PSD_gatherer_fg_burst,2);
    avg_PSD_gatherer_fg_burst=movmean(avg_PSD_gatherer_fg_burst,avg_window);
    plot(freqVals,avg_PSD_gatherer_fg_burst,'LineWidth',2,'Color',color_orange);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    random_avg_PSD_gatherer_sg_burst=mean(random_PSD_gatherer_sg_burst,2);
    random_avg_PSD_gatherer_sg_burst=movmean(random_avg_PSD_gatherer_sg_burst,avg_window);
    % dotted blue line for random slow gamma burst
    plot(freqVals,random_avg_PSD_gatherer_sg_burst,'--b','LineWidth',1.5);
    random_avg_PSD_gatherer_fg_burst=mean(random_PSD_gatherer_fg_burst,2);
    random_avg_PSD_gatherer_fg_burst=movmean(random_avg_PSD_gatherer_fg_burst,avg_window);
    % dotted orange line for random fast gamma burst
    plot(freqVals,random_avg_PSD_gatherer_fg_burst,'--','Color',color_orange,'LineWidth',1.5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    % Draw vertical lines to indicate slow and fast gamma frequency ranges
    xline(slow_gamma_freq(1),'--b','LineWidth',1.5);
    xline(slow_gamma_freq(2),'--b','LineWidth',1.5);
    xline(fast_gamma_freq(1),'--','Color',color_orange,'LineWidth',1.5);
    xline(fast_gamma_freq(2),'--','Color',color_orange,'LineWidth',1.5);
    legend('Slow \gamma','Fast \gamma','Slow \gamma (Random)','Fast \gamma (Random)');
    xlim([10 70]);
    ylim([-40 20]);
end