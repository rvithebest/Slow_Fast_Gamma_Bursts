clc;clear
f=figure;
f.WindowState="Maximized";
plotHandles_1=getPlotHandles(2,1,[0.01 0.08 0.3 0.9],0.07,0.09,0);
plotHandles_2=getPlotHandles(2,4,[0.4 0.08 0.58 0.9],0.05,0.09,0);
load('Model_sim_SG_theta_64_4_input_14\MP_analysis_results.mat');
LFP_sg=lfp_decimated;
load('Model_sim_FG_theta_16_1_input_7\MP_analysis_results.mat');
% Representative parametric conditions
% For confirming the generation of sg and fg oscillations
LFP_fg=lfp_decimated;
timeVals=timeVals_decimated;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_1(1,1));
img = imread('model_schematic_1.jpg');
imagesc(img);
axis image off;
axis off;
subplot(plotHandles_1(2,1));
img = imread('model_schematic_2.jpg');
imagesc(img);
axis image off;
axis off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sg_freq=[20 40]; fg_freq=[40 80];
Model_TF_plot(LFP_sg,timeVals,1,plotHandles_2,sg_freq);
Model_TF_plot(LFP_fg,timeVals,2,plotHandles_2,fg_freq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_E_sg=16; theta_I_sg=1; input_sg=14;
theta_E_fg=16; theta_I_fg=1; input_fg=7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sg_idx=find_idx(theta_E_sg, theta_I_sg, input_sg);
fg_idx=find_idx(theta_E_fg, theta_I_fg, input_fg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_2(1,2));
% PSD of LFP_sg and LFP_fg
Fs=1/(timeVals(2)-timeVals(1));
stimPeriod=[0.25 1.25];
time_idx=find(timeVals>=stimPeriod(1) & timeVals<=stimPeriod(2));
TW=2;
K=3;
params.Fs=Fs;
params.tapers=[TW K];
params.fpass=[0 250];
params.trialave=1;
params.pad=-1; % Disable zero padding
params.err=[2 0.05];
[S_sg,f_sg]=mtspectrumc(LFP_sg(:,time_idx)',params);
[S_fg,f_fg]=mtspectrumc(LFP_fg(:,time_idx)',params);
plot(f_sg,log10(S_sg),'LineWidth',2,'Color','b');
hold on;
color_orange=[1 0.5 0]; % Define the orange color
color_violet = [0.54 0.17 0.89]; % RGB triplet for violet
plot(f_fg,log10(S_fg),'LineWidth',2,'Color',color_orange);
xline([sg_freq(1) sg_freq(1)],'Color','b','LineStyle','--','LineWidth',2);
xline([sg_freq(2) sg_freq(2)],'Color',color_violet,'LineStyle','--','LineWidth',2);
xline([fg_freq(1) fg_freq(1)],'Color',color_violet,'LineStyle','--','LineWidth',2);
xline([fg_freq(2) fg_freq(2)],'Color',color_orange,'LineStyle','--','LineWidth',2);
xlim([0 100]);
xlabel("Frequency (Hz)");
ylabel("Raw Power (log)")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MP_plot_model_sim(plotHandles_2,sg_idx,fg_idx,sg_freq,fg_freq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set_axis_ticks_fontsize(plotHandles_2,22,16,1);
% set_axis_ticks_fontsize(plotHandles_2,22,16,2);
set_axis_ticks_fontsize(plotHandles_2,16,14,1);
set_axis_ticks_fontsize(plotHandles_2,16,14,2);
labels = {'A','B','C','D','E','F','G','H','I',"J"};
x_positions = [0.01, 0.35, 0.51, 0.68, 0.83];
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
function idx = find_idx(theta_E, theta_I, input)
    % Find the index based on the parameters
    theta_E_vals=[1,4,16,64];
    theta_I_vals=[1,4,16,64];
    input_vals=[7,8,9,10,11,12,13,14];
    % 1,1,7- idx 1
    % 1,1,8- idx 2
    % 1,1,15- idx 9
    % 1,4,7- idx 10
    idx_E=find(theta_E_vals==theta_E);
    idx_I=find(theta_I_vals==theta_I);
    idx_input=find(input_vals==input);
    if isempty(idx_E) || isempty(idx_I) || isempty(idx_input)
        error('Invalid parameters provided.');
    end
    idx = ((idx_E-1)*(length(theta_I_vals)*length(input_vals))) + ((idx_I-1)*length(input_vals)) + idx_input;
end