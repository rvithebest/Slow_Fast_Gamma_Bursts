clc;clear
f=figure;
f.WindowState="Maximized";
plotHandles_1=getPlotHandles(2,1,[0.08 0.08 0.8 0.8],0.15,0.15,0);
slow_gamma_freq=[20 40]; fast_gamma_freq=[40 80];
theta_E_vals=[1,4,16,64];
theta_I_vals=[1,4,16,64];
input_vals=[7,8,9,10,11,12,13,14];
length_vals_sg=cell(4,4);
length_vals_fg=cell(4,4);
onset_vals_fg=cell(4,4);
onset_vals_sg=cell(4,4);
idx=1;
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
            %%%%%%%%%%%%% Burst Length %%%%%%%%%%%%%%%
            % [Power_fg] = get_Power_model(lfp_decimated, timeVals_decimated, stimulusPeriodS, 2);
            [Thresh_vals_fg,Width_vals_fg] = get_Threshold_model(lfp_decimated, timeVals_decimated, stimulusPeriodS,2);
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
                onset_time_temp_trial_fg(onset_time_temp_trial_fg<0)=0; % Ensure onset times are non-negative (after stimulus presentation)
                onset_idx=(find((onset_time_temp_trial_fg)==min(onset_time_temp_trial_fg)));
                onset_temp_all_trials_fg=[onset_temp_all_trials_fg,onset_time_temp_trial_fg(onset_idx)];
            end
            fg_l=length_temp_all_trials_fg;
            if ~isempty(fg_l)
                length_vals_fg{i,j}=[length_vals_fg{i,j},fg_l];
            end
            onset_vals_fg{i,j}=[onset_vals_fg{i,j},onset_temp_all_trials_fg];
            %%%%%%%%%%%%%%%%% SG %%%%%%%%%%%%%%%%%%%%
            folder_name="Model_sim_SG_theta_"+string(theta_E)+"_"+string(theta_I)+"_input_"+string(input);
            % load simulationsResults.mat inside the folder
            load(folder_name+"/MP_analysis_results.mat");
            gamma_freq=slow_gamma_freq;
            %%%%%%%%%%%%%%%% Burst Length %%%%%%%%%%%%%%%
            % [Power_sg]= get_Power_model(lfp_decimated, timeVals_decimated, stimulusPeriodS,1);
            [Thresh_vals_sg, Width_vals_sg] = get_Threshold_model(lfp_decimated, timeVals_decimated, stimulusPeriodS,1);
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
                onset_time_temp_trial_sg(onset_time_temp_trial_sg<0)=0; % Ensure onset times are non-negative (after stimulus presentation)
                onset_idx=(find((onset_time_temp_trial_sg)==min(onset_time_temp_trial_sg)));
                onset_temp_all_trials_sg=[onset_temp_all_trials_sg,onset_time_temp_trial_sg(onset_idx)];
            end
            sg_l=length_temp_all_trials_sg;
            if ~isempty(sg_l)
               length_vals_sg{i,j}=[length_vals_sg{i,j},sg_l];
            end
            onset_vals_sg{i,j}=[onset_vals_sg{i,j},onset_temp_all_trials_sg];
            idx=idx+1;
        end
    end
end
median_length_vals_sg=zeros(4,4);
median_length_vals_fg=zeros(4,4);
mean_onset_vals_sg=zeros(4,4);
mean_onset_vals_fg=zeros(4,4);
for i=1:4
for j=1:4
    median_length_vals_sg(i,j)=median(length_vals_sg{i,j});
    median_length_vals_fg(i,j)=median(length_vals_fg{i,j});
    % mean_onset_vals_sg(i,j)=mean(onset_vals_sg{i,j});
    % mean_onset_vals_fg(i,j)=mean(onset_vals_fg{i,j});
end
end
% Plot heatmaps of median burst lengths for slow gamma and fast gamma
subplot(plotHandles_1(1,1));
imagesc(median_length_vals_sg);
% revert the y-axis
set(gca,'YDir','normal');
c=colorbar;
% set climits
caxis([0.2 1.4]);
% Y ticks- theta_I (1,4,16,64)
theta_I_axis=categorical({'1','4','16','64'});
set(gca,'XTick',1:4,'YTickLabel',theta_I_axis);
theta_E_axis=categorical({'1','4','16','64'});
set(gca,'YTick',1:4,'XTickLabel',theta_E_axis);
title('Median Burst Length - Slow Gamma');
xlabel('\theta_I');
ylabel('\theta_E');
subplot(plotHandles_1(2,1));
imagesc(median_length_vals_fg);
% revert the y-axis
set(gca,'YDir','normal');
c=colorbar;
% set climits
caxis([0.2 1.4]);
% Y ticks- theta_I (1,4,16,64)
theta_I_axis=categorical({'1','4','16','64'});
set(gca,'XTick',1:4,'YTickLabel',theta_I_axis);
theta_E_axis=categorical({'1','4','16','64'});
set(gca,'YTick',1:4,'XTickLabel',theta_E_axis);
title('Median Burst Length - Fast Gamma');
xlabel('\theta_I');
ylabel('\theta_E');
set_axis_ticks_fontsize(plotHandles_1,16,14,1);
set_axis_ticks_fontsize(plotHandles_1,16,14,2);
labels = {'A','B'};
x_positions = [0.01];
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
% subplot(2,2,3);
% imagesc(mean_onset_vals_sg);
% % revert the y-axis 
% set(gca,'YDir','normal');
% colorbar;
% % Y ticks- theta_I (1,4,16,64)
% theta_I_axis=categorical({'1','4','16','64'});
% set(gca,'XTick',1:4,'YTickLabel',theta_I_axis);
% theta_E_axis=categorical({'1','4','16','64'});
% set(gca,'YTick',1:4,'XTickLabel',theta_E_axis);
% title('Mean Onset Time - Slow Gamma');
% xlabel('\theta_I');
% ylabel('\theta_E');
% subplot(2,2,4);
% imagesc(mean_onset_vals_fg);
% % revert the y-axis
% set(gca,'YDir','normal');
% colorbar;
% % Y ticks- theta_I (1,4,16,64)
% theta_I_axis=categorical({'1','4','16','64'});
% set(gca,'XTick',1:4,'YTickLabel',theta_I_axis);
% theta_E_axis=categorical({'1','4','16','64'});
% set(gca,'YTick',1:4,'XTickLabel',theta_E_axis);
% title('Mean Onset Time - Fast Gamma');
