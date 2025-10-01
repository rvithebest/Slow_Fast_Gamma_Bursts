function MP_plot_model_sim(plotHandles,sg_idx,fg_idx,slow_gamma_freq,fast_gamma_freq)
    theta_E_vals=[1,4,16,64];
    theta_I_vals=[1,4,16,64];
    input_vals=[7,8,9,10,11,12,13,14];
    N=(length(theta_E_vals)*length(theta_I_vals)*length(input_vals));
    median_sg_vals=zeros(1,N);
    median_fg_vals=zeros(1,N);
    mean_o_sg_vals=zeros(1,N);
    mean_o_fg_vals=zeros(1,N);
    power_sg_vals=zeros(1,N);
    power_fg_vals=zeros(1,N);
    length_all_elec_sg=[];
    length_all_elec_fg=[];
    idx=1;
    counter=0;
    for i=1:length(theta_E_vals)
        for j=1:length(theta_I_vals)
            for k=1:length(input_vals)
                theta_E=theta_E_vals(i);
                theta_I=theta_I_vals(j);
                input=input_vals(k);
                %%%%Z%%%%%%%%%% FG %%%%%%%%%%%%%%%%%%%%
                folder_name="Model_sim_FG_theta_"+string(theta_E)+"_"+string(theta_I)+"_input_"+string(input);
                % load simulationsResults.mat inside the folder
                load(folder_name+"/MP_analysis_results.mat");
                gamma_freq=fast_gamma_freq;
                [Power_fg] = get_Power_model(lfp_decimated, timeVals_decimated, stimulusPeriodS, 2);
                [Thresh_vals_fg,Width_vals_fg] = get_Threshold_model(lfp_decimated, timeVals_decimated, stimulusPeriodS,2);
                power_fg_vals(idx)=Power_fg;
                %%%%%%%%%%%%% Burst Length %%%%%%%%%%%%%%%
                [length_temp_fg,~,time_center_temp_fg,gabor_temp,header_temp,~]= getBurstLengthMP_model(lfp_decimated,timeVals_decimated,displayFlag,stimulusPeriodS,baselinePeriodS,gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp,Thresh_vals_fg,Width_vals_fg);
                length_temp_all_trials_fg=[];
                onset_temp_all_trials_fg=[];
                for ii=1:length(length_temp_fg)
                    if isempty(length_temp_fg{ii}')
                        continue;
                    end
                    reject_idx=find((length_temp_fg{ii}')>1.8);
                    length_temp_fg{ii}(reject_idx)=[];
                    time_center_temp_fg{ii}(reject_idx)=[];
                    length_temp_all_trials_fg=[length_temp_all_trials_fg,length_temp_fg{ii}'];
                    onset_time_temp_trial_fg=time_center_temp_fg{ii}'-((length_temp_fg{ii}')*0.5);
                    onset_time_temp_trial_fg(onset_time_temp_trial_fg<0)=[]; % Ensure onset times are non-negative (after stimulus presentation)
                    onset_idx=(find((onset_time_temp_trial_fg)==min(onset_time_temp_trial_fg)));
                    onset_temp_all_trials_fg=[onset_temp_all_trials_fg,onset_time_temp_trial_fg(onset_idx)];
                end
                fg_l=length_temp_all_trials_fg;
                if (((~isempty(fg_l)))&&(~isempty(onset_temp_all_trials_fg)))
                    median_fg=median(fg_l);
                    median_fg_vals(idx)=median_fg;
                    mean_o_fg_vals(idx)=mean(onset_temp_all_trials_fg);
                else
                    mean_o_fg_vals(idx)=0;
                    median_fg_vals(idx)=0;
                    power_fg_vals(idx)=0;
                    median_fg=0;
                end
                %%%%%%%%%%%%%%%%% SG %%%%%%%%%%%%%%%%%%%%
                folder_name="Model_sim_SG_theta_"+string(theta_E)+"_"+string(theta_I)+"_input_"+string(input);
                % load simulationsResults.mat inside the folder
                load(folder_name+"/MP_analysis_results.mat");
                gamma_freq=slow_gamma_freq;
                [Power_sg]= get_Power_model(lfp_decimated, timeVals_decimated, stimulusPeriodS, 1);
                [Thresh_vals_sg,Width_vals_sg] = get_Threshold_model(lfp_decimated, timeVals_decimated, stimulusPeriodS,1);
                power_sg_vals(idx)=Power_sg;
                %%%%%%%%%%%%%%%% Burst Length %%%%%%%%%%%%%%%
                [length_temp_sg,~,time_center_temp_sg,gabor_temp,header_temp,~]= getBurstLengthMP_model(lfp_decimated,timeVals_decimated,displayFlag,stimulusPeriodS,baselinePeriodS,gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp,Thresh_vals_sg,Width_vals_sg);
                length_temp_all_trials_sg=[];
                onset_temp_all_trials_sg=[];
                for ii=1:length(length_temp_sg)
                    if isempty(length_temp_sg{ii}')
                        continue;
                    end
                    reject_idx=find((length_temp_sg{ii}')>1.8);
                    length_temp_sg{ii}(reject_idx)=[];
                    time_center_temp_sg{ii}(reject_idx)=[];
                    length_temp_all_trials_sg=[length_temp_all_trials_sg,length_temp_sg{ii}'];
                    onset_time_temp_trial_sg=time_center_temp_sg{ii}'-((length_temp_sg{ii}')*0.5);
                    onset_time_temp_trial_sg(onset_time_temp_trial_sg<0)=[]; % Ensure onset times are non-negative (after stimulus presentation)
                    onset_idx=(find((onset_time_temp_trial_sg)==min(onset_time_temp_trial_sg)));
                    onset_temp_all_trials_sg=[onset_temp_all_trials_sg,onset_time_temp_trial_sg(onset_idx)];
                end
                sg_l=length_temp_all_trials_sg;
                if (((~isempty(sg_l)))&&(~isempty(onset_temp_all_trials_sg)))
                   median_sg=median(sg_l);
                   median_sg_vals(idx)=median_sg;
                   mean_o_sg_vals(idx)=mean(onset_temp_all_trials_sg);
                else 
                    median_sg_vals(idx)=0;
                    mean_o_sg_vals(idx)=0;
                    power_sg_vals(idx)=0;
                    median_fg_vals(idx)=0;
                    mean_o_fg_vals(idx)=0;
                    power_fg_vals(idx)=0;
                    median_sg=0;
                    median_fg=0;
                end
                if (median_fg_vals(idx)==0)
                    power_sg_vals(idx)=0;
                    median_sg_vals(idx)=0;
                    mean_o_sg_vals(idx)=0;
                end
                if median_sg < median_fg
                    counter=counter+1;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if idx== sg_idx
                    sgl_plot=[];
                    sgl_plot=sg_l;
                end
                if idx== fg_idx
                    fgl_plot=[];
                    fgl_plot=fg_l;
                end
                % If sg_l or fg_l has less than 40 bursts, we do not analyze it
                % if length(sg_l)<40 || length(fg_l)<40
                %     median_sg_vals(idx)=0;
                %     median_fg_vals(idx)=0;
                %     power_sg_vals(idx)=0;
                %     power_fg_vals(idx)=0;
                % end
                length_all_elec_sg=[length_all_elec_sg,sg_l];
                length_all_elec_fg=[length_all_elec_fg,fg_l];
                idx=idx+1;
            end
        end
    end
    disp(counter);
    % Reject all the null values
    median_sg_vals(median_sg_vals==0)=[];
    median_fg_vals(median_fg_vals==0)=[];
    power_sg_vals(power_sg_vals==0)=[];
    power_fg_vals(power_fg_vals==0)=[];
    mean_o_sg_vals(mean_o_sg_vals==0)=[];
    mean_o_fg_vals(mean_o_fg_vals==0)=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(2,2));
    hold on;
    % xline(median(fgl_plot),'--','Color',[1 0.5 0],'LineWidth',2);
    % xline(median(sgl_plot),'--','Color','b','LineWidth',2);
    % histogram_all_burst(sgl_plot,fgl_plot,2,10);
    histogram_all_burst(length_all_elec_sg,length_all_elec_fg,2,10);
    xlim([0 2]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subplot(plotHandles(2,2));
    % % Only for equal inputs- (8,8) for both the models
    % plot_idx=2:length(input_vals):N;
    % violin_swarm_plot_paired(median_sg_vals(plot_idx),median_fg_vals(plot_idx),0);
    subplot(plotHandles(1,3));
    violin_swarm_plot_paired(mean_o_sg_vals, mean_o_fg_vals,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(2,3));
    violin_swarm_plot_paired(median_sg_vals,median_fg_vals,0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(1,4));
    [matched_sg_indices,matched_fg_indices]=power_matching_hist(power_sg_vals,power_fg_vals);
    scatter_plot_power_burst(power_sg_vals,power_fg_vals,median_sg_vals ...
        ,median_fg_vals,matched_sg_indices,matched_fg_indices,2);
    hold on;
    % plot(power_sg_vals(sg_idx),median_sg_vals(sg_idx),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);
    % plot(power_fg_vals(fg_idx),median_fg_vals(fg_idx),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(plotHandles(2,4));
    matched_sg_lengths= median_sg_vals(matched_sg_indices);
    matched_fg_lengths= median_fg_vals(matched_fg_indices);
    violin_swarm_plot(matched_sg_lengths, matched_fg_lengths);
    save("Model_analysis_data_v2.mat","length_all_elec_fg","length_all_elec_sg","mean_o_fg_vals","mean_o_sg_vals","matched_fg_lengths","matched_sg_lengths","matched_fg_indices","matched_sg_indices","power_fg_vals","power_sg_vals","median_fg_vals","median_sg_vals");
end