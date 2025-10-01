function plot_conn(method, plotHandles_1,plotHandles_2,idx)
    for Monkey_num=1:2
        clearvars -except Monkey_num f plotHandles_1 plotHandles_2 method idx
        parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Monkey_num==1
            % Monkey- alpaH
            load((fullfile(parent_file_path,'alpaH_info','parameterCombinations.mat')))
            load(fullfile(parent_file_path,'alpaH_info','badTrials.mat'));
            load(fullfile(parent_file_path,'alpaH_info','alpaHMicroelectrodeRFData.mat'));
            LFP_data_file=dir(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H'));
            LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'}));
            LFP_data_file = natsortfiles({LFP_data_file.name});
            % slow_gamma_freq=[20 34];
            % fast_gamma_freq=[34 65];
            if strcmp(method,'PLV')
                load('M1_conn_analysis_PLV.mat');
                slow_gamma_freq=[28 28];
                fast_gamma_freq=[42 42];
            elseif strcmp(method,'WPLI')
                load('M1_conn_analysis_WPLI.mat');
                slow_gamma_freq=[28 28];
                fast_gamma_freq=[46 46];
            else
                error('Invalid method specified. Use "PLV" or "WPLI".');
            end
            selected_elec=highRMSElectrodes;
            selected_elec_LFP=selected_elec(1:77);
        elseif Monkey_num==2
            % Monkey- kesariH
            load((fullfile(parent_file_path,'kesariH_info','parameterCombinations.mat')))
            load(fullfile(parent_file_path,'kesariH_info','badTrials_kesari.mat'));
            load(fullfile(parent_file_path,'kesariH_info','kesariHMicroelectrodeRFData_Two.mat'));
            LFP_data_file=dir(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H'));
            LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'}));
            LFP_data_file = natsortfiles({LFP_data_file.name});
            % slow_gamma_freq=[20 40];
            % fast_gamma_freq=[40 65];
            if strcmp(method,'PLV')
                load('M2_conn_analysis_PLV.mat');
                slow_gamma_freq=[20 20];
                fast_gamma_freq=[42 42];
            elseif strcmp(method,'WPLI')
                load('M2_conn_analysis_WPLI.mat');
                slow_gamma_freq=[26 26];
                fast_gamma_freq=[50 50];
            else
                error('Invalid method specified. Use "PLV" or "WPLI".');
            end
            selected_elec=highRMSElectrodes;
            selected_elec_LFP=selected_elec(1:31);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dist=[400,1200,3200,4800];
        load('timeVals.mat');
        num_elec=length(selected_elec_LFP);
        % conn_stim is 77x77x51 double
        freq_axis=freq_stim.freq;
        % freq_axis=freq_vals_ST{1,1};
        % Find idx where frequency is in slow gamma range
        slow_gamma_idx=find(freq_axis>=slow_gamma_freq(1) & freq_axis<=slow_gamma_freq(2));
        % Find idx where frequency is in fast gamma range
        fast_gamma_idx=find(freq_axis>=fast_gamma_freq(1) & freq_axis<=fast_gamma_freq(2));
        dist_group_1_freq=[]; % Frequency group 1
        dist_group_2_freq=[]; % Frequency group 2
        dist_group_3_freq=[]; % Frequency group 3
        dist_group_4_freq=[]; % Frequency group 4
        dist_group_5_freq=[]; % Frequency group 5
        dist_group_1_sg=[]; 
        dist_group_2_sg=[];
        dist_group_3_sg=[]; 
        dist_group_4_sg=[]; 
        dist_group_1_fg=[]; 
        dist_group_2_fg=[]; 
        dist_group_3_fg=[]; 
        dist_group_4_fg=[];
        dist_group_5_sg=[]; 
        dist_group_5_fg=[];
        % N Choose 2 combinations of electrodes- (N*(N-1))/2
        % Loop through each unique pair of electrodes
        for i=1:(num_elec-1)
           distance_array=getElectrodeDistance(selected_elec_LFP,selected_elec_LFP(i));
           for j=i+1:(num_elec)
            current_distance=distance_array(j);
            conn_temp_stim=squeeze(conn_stim(i,j,:));
            conn_temp_bl=squeeze(conn_bl(i,j,:));
            conn_temp_sg=max(conn_temp_stim(slow_gamma_idx));
            conn_temp_fg=max(conn_temp_stim(fast_gamma_idx));
            if current_distance>0 && current_distance<=dist(1)
                dist_group_1_freq=[dist_group_1_freq;conn_temp_stim'];
                dist_group_1_sg=[dist_group_1_sg;conn_temp_sg];
                dist_group_1_fg=[dist_group_1_fg;conn_temp_fg];
            elseif current_distance>dist(1) && current_distance<=dist(2)
                dist_group_2_freq=[dist_group_2_freq;conn_temp_stim'];
                dist_group_2_sg=[dist_group_2_sg;conn_temp_sg];
                dist_group_2_fg=[dist_group_2_fg;conn_temp_fg];
            elseif current_distance>dist(2) && current_distance<=dist(3)
                dist_group_3_freq=[dist_group_3_freq;conn_temp_stim'];
                dist_group_3_sg=[dist_group_3_sg;conn_temp_sg];
                dist_group_3_fg=[dist_group_3_fg;conn_temp_fg];
            elseif current_distance>dist(3) && current_distance<=dist(4)
                dist_group_4_freq=[dist_group_4_freq;conn_temp_stim'];
                dist_group_4_sg=[dist_group_4_sg;conn_temp_sg];
                dist_group_4_fg=[dist_group_4_fg;conn_temp_fg];
            else
                dist_group_5_freq=[dist_group_5_freq;conn_temp_stim'];
                dist_group_5_sg=[dist_group_5_sg;conn_temp_sg];
                dist_group_5_fg=[dist_group_5_fg;conn_temp_fg];
            end
           end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(plotHandles_1(Monkey_num,idx));
        % red
        color_group_1= 'r';
        % dark green
        color_group_2= [0 0.5 0];
        % magenta
        color_group_3= 'm';
        % purple
        color_group_4= [0.5 0 0.5];
        % Mean analysis
        plot(freq_axis,mean(dist_group_1_freq),'LineWidth',2,'Color',color_group_1);
        hold on;
        plot(freq_axis,mean(dist_group_2_freq),'LineWidth',2,'Color',color_group_2);    
        plot(freq_axis,mean(dist_group_3_freq),'LineWidth',2,'Color',color_group_3);
        plot(freq_axis,mean(dist_group_4_freq),'LineWidth',2,'Color',color_group_4);
        xline(slow_gamma_freq(1),'--','Color','b','LineWidth',2);
        xline(fast_gamma_freq(1),'--','Color',[1 0.5 0],'LineWidth',2);
        xlim([10 80]);
        if Monkey_num==1
           hold on;
           if strcmp(method,'PLV')
                % Text positions
                x_text = 55;
                y_vals = [0.67, 0.75, 0.83, 0.91];
                dist_1_string=string(dist(1));
                dist_2_string="("+string(dist(1))+","+string(dist(2))+"]";
                dist_3_string="("+string(dist(2))+","+string(dist(3))+"]";
                dist_4_string="("+string(dist(3))+","+string(dist(4))+"]";
                % Plot color-coded texts
                text(x_text, y_vals(4), dist_1_string, 'Color', color_group_1, 'FontSize', 18);
                text(x_text, y_vals(3), dist_2_string, 'Color', color_group_2, 'FontSize', 18); 
                text(x_text, y_vals(2), dist_3_string, 'Color', color_group_3, 'FontSize', 18);
                text(x_text, y_vals(1), dist_4_string, 'Color', color_group_4, 'FontSize', 18);
                
                % Compute bounding box dimensions
                margin_x = 2;  % horizontal padding
                margin_y = 0.05; % vertical padding
                
                % Rectangle lower-left corner
                x_box = x_text - margin_x;
                y_box = y_vals(1) - margin_y;
                
                % Rectangle height and width
                box_height = (y_vals(end) - y_vals(1)) + 2*margin_y;
                box_width  = 25;  % adjust as needed
                
                % Draw box
                rectangle('Position', [x_box, y_box, box_width, box_height], ...
                          'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
                % Optionally set x-ticks to empty to rely only on annotations
                set(gca, 'XTickLabel', []);
            end
            if strcmp(method,'WPLI')
                % Text positions
                x_text = 70;
                y_vals = [0.59,0.64, 0.69, 0.74];
                N1= size(dist_group_1_freq, 1);
                N2= size(dist_group_2_freq, 1);
                N3= size(dist_group_3_freq, 1);
                N4= size(dist_group_4_freq, 1);
                string_1='N='+string(N1);
                % Plot color-coded texts
                text(x_text, y_vals(4), string_1, 'Color', color_group_1, 'FontSize', 14);
                text(x_text, y_vals(3), string(N2), 'Color', color_group_2, 'FontSize', 14);
                text(x_text, y_vals(2), string(N3), 'Color', color_group_3, 'FontSize', 14);
                text(x_text, y_vals(1), string(N4), 'Color', color_group_4, 'FontSize', 14);
                % Compute bounding box dimensions
                margin_x = 2;  % horizontal padding
                margin_y = 0.05; % vertical padding
                
                % Rectangle lower-left corner
                x_box = x_text - margin_x;
                y_box = y_vals(1) - margin_y;
                
                % Rectangle height and width
                box_height = (y_vals(end) - y_vals(1)) + 2*margin_y;
                box_width  = 30;  % adjust as needed
                
                % Draw box
                rectangle('Position', [x_box, y_box, box_width, box_height], ...
                          'EdgeColor', 'w', 'LineWidth', 1.5, 'LineStyle', '-');
                % Optionally set x-ticks to empty to rely only on annotations
                set(gca, 'XTickLabel', []);
            end
    
        else 
           xlabel('Frequency (Hz)');
            if strcmp(method,'WPLI')
                % Text positions
                x_text = 70;
                y_vals = [0.54,0.62, 0.7, 0.78];
                N1= size(dist_group_1_freq, 1);
                N2= size(dist_group_2_freq, 1);
                N3= size(dist_group_3_freq, 1);
                N4= size(dist_group_4_freq, 1);
                string_1='N='+string(N1);
                % Plot color-coded texts
                text(x_text, y_vals(4), string_1, 'Color', color_group_1, 'FontSize', 14);
                text(x_text, y_vals(3), string(N2), 'Color', color_group_2, 'FontSize', 14);
                text(x_text, y_vals(2), string(N3), 'Color', color_group_3, 'FontSize', 14);
                text(x_text, y_vals(1), string(N4), 'Color', color_group_4, 'FontSize', 14);
                % Compute bounding box dimensions
                margin_x = 2;  % horizontal padding
                margin_y = 0.05; % vertical padding
                
                % Rectangle lower-left corner
                x_box = x_text - margin_x;
                y_box = y_vals(1) - margin_y;
                
                % Rectangle height and width
                box_height = (y_vals(end) - y_vals(1)) + 2*margin_y;
                box_width  = 30;  % adjust as needed
                
                % Draw box
                rectangle('Position', [x_box, y_box, box_width, box_height], ...
                          'EdgeColor', 'w', 'LineWidth', 1.5, 'LineStyle', '-');
            end
        end
        if strcmp(method,'PLV')
            ylabel("PLV");
        elseif strcmp(method,'WPLI')
            ylabel("WPLI");
        end
        % Summary plots for slow and fast gamma
        subplot(plotHandles_2(Monkey_num,idx));
        avg_dist_group_1_fg=mean(dist_group_1_fg);
        avg_dist_group_2_fg=mean(dist_group_2_fg);
        avg_dist_group_3_fg=mean(dist_group_3_fg);
        avg_dist_group_4_fg=mean(dist_group_4_fg);
        avg_dist_group_1_sg=mean(dist_group_1_sg);
        avg_dist_group_2_sg=mean(dist_group_2_sg);
        avg_dist_group_3_sg=mean(dist_group_3_sg);
        avg_dist_group_4_sg=mean(dist_group_4_sg);
        SEM_dist_group_1_fg=std(dist_group_1_fg)/sqrt(length(dist_group_1_fg));
        SEM_dist_group_2_fg=std(dist_group_2_fg)/sqrt(length(dist_group_2_fg));
        SEM_dist_group_3_fg=std(dist_group_3_fg)/sqrt(length(dist_group_3_fg));
        SEM_dist_group_4_fg=std(dist_group_4_fg)/sqrt(length(dist_group_4_fg));
        SEM_dist_group_1_sg=std(dist_group_1_sg)/sqrt(length(dist_group_1_sg));
        SEM_dist_group_2_sg=std(dist_group_2_sg)/sqrt(length(dist_group_2_sg));
        SEM_dist_group_3_sg=std(dist_group_3_sg)/sqrt(length(dist_group_3_sg));
        SEM_dist_group_4_sg=std(dist_group_4_sg)/sqrt(length(dist_group_4_sg));
        errorbar([1 2 3 4],[avg_dist_group_1_sg avg_dist_group_2_sg avg_dist_group_3_sg avg_dist_group_4_sg],...
         [SEM_dist_group_1_sg SEM_dist_group_2_sg SEM_dist_group_3_sg SEM_dist_group_4_sg],'o-','LineWidth',2,'MarkerSize',10,'Color','b');
        color_orange=[1 0.5 0];
        hold on;
        errorbar([1 2 3 4],[avg_dist_group_1_fg avg_dist_group_2_fg avg_dist_group_3_fg avg_dist_group_4_fg],...
             [SEM_dist_group_1_fg SEM_dist_group_2_fg SEM_dist_group_3_fg SEM_dist_group_4_fg],'o-','LineWidth',2,'MarkerSize',10,'Color',color_orange);
        hold on;
        plot([1],avg_dist_group_1_sg,'o','MarkerSize',10,'LineWidth',2,'Color',color_group_1,'MarkerFaceColor',color_group_1);
        plot([2],avg_dist_group_2_sg,'o','MarkerSize',10,'LineWidth',2,'Color',color_group_2,'MarkerFaceColor',color_group_2);
        plot([3],avg_dist_group_3_sg,'o','MarkerSize',10,'LineWidth',2,'Color',color_group_3,'MarkerFaceColor',color_group_3);
        plot([4],avg_dist_group_4_sg,'o','MarkerSize',10,'LineWidth',2,'Color',color_group_4,'MarkerFaceColor',color_group_4);
        plot([1],avg_dist_group_1_fg,'o','MarkerSize',10,'LineWidth',2,'Color',color_group_1,'MarkerFaceColor',color_group_1);
        plot([2],avg_dist_group_2_fg,'o','MarkerSize',10,'LineWidth',2,'Color',color_group_2,'MarkerFaceColor',color_group_2);
        plot([3],avg_dist_group_3_fg,'o','MarkerSize',10,'LineWidth',2,'Color',color_group_3,'MarkerFaceColor',color_group_3);
        plot([4],avg_dist_group_4_fg,'o','MarkerSize',10,'LineWidth',2,'Color',color_group_4,'MarkerFaceColor',color_group_4);
        % 2 way ANOVA analysis
        allData=[dist_group_1_sg;dist_group_2_sg;dist_group_3_sg;dist_group_4_sg;...
            dist_group_1_fg;dist_group_2_fg;dist_group_3_fg;dist_group_4_fg];
        group_labels=[repmat("1",length(dist_group_1_sg),1);repmat("2",length(dist_group_2_sg),1);...
            repmat("3",length(dist_group_3_sg),1);repmat("4",length(dist_group_4_sg),1);...
            repmat("1",length(dist_group_1_fg),1);repmat("2",length(dist_group_2_fg),1);...
            repmat("3",length(dist_group_3_fg),1);repmat("4",length(dist_group_4_fg),1)];
        cond_labels=[repmat("Slow \gamma",length(dist_group_1_sg)+length(dist_group_2_sg)+...
            length(dist_group_3_sg)+length(dist_group_4_sg),1);...
            repmat("Fast \gamma",length(dist_group_1_fg)+length(dist_group_2_fg)+...
            length(dist_group_3_fg)+length(dist_group_4_fg),1)];
        [p,tble,stats] = anovan(allData, {group_labels, cond_labels},'model','interaction','varnames',{'Distance group','Rhythm type'});
        df_Int = tble{4,3}; F_Int = tble{4,6}; p_Int = tble{4,7}; df_error  = tble{end,3};
        % Legend- for the colors (customized)- interelec distances
        % dist(1), (dist(1),dist(2)], (dist(2),dist(3)], (dist(3),dist(4)] (font color should be same as the color of the group)
        % Add custom color-coded text annotations near x-ticks
        % Optionally set x-ticks to empty to rely only on annotations
        set(gca, 'XTickLabel', []);
        if Monkey_num==2
           xlabel("Distance group");
        else 
            dummy1 = plot(NaN, NaN, '-', 'LineWidth', 2, 'Color', 'b');
            dummy2 = plot(NaN, NaN, '-', 'LineWidth', 2, 'Color', color_orange);
            lgd = legend([dummy1, dummy2], ...
            {'Slow \gamma','Fast \gamma'}, ...
            'Location', 'northeast'); 
        end
        if strcmp(method,'PLV')
            ylabel("Peak PLV");
        elseif strcmp(method,'WPLI')
            ylabel("Peak WPLI");
        end
    end
end