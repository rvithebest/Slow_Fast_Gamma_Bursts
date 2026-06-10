clc;clear;close all;
length_accumulator_omp_gear=cell(8,1);
length_accumulator_omp_gear_tf_l=cell(8,1);
length_accumulator_omp_gear_tf_h=cell(8,1);
gabor_accmulator_omp_gear=cell(8,1);
length_accumulator_MP=cell(8,1);
length_accumulator_MP_tf_l=cell(8,1);
length_accumulator_MP_tf_h=cell(8,1);
gabor_accmulator_MP=cell(8,1);
length_accumulator_hilbert=cell(8,1);
length_accumulator_hilbert_tf_l=cell(8,1);
length_accumulator_hilbert_tf_h=cell(8,1);
length_accumulator_wavelet=cell(8,1);
length_accumulator_wavelet_tf_l=cell(8,1);
length_accumulator_wavelet_tf_h=cell(8,1);
analogData_accumulator=cell(8,1);
analogData0_accumulator=cell(8,1); % Spontaneous signal before burst injection
timeVals_accumulator=cell(8,1); % Time values for the synthetic data
numBurst_accumulator=cell(8,1); 
freqBurst_accumulator=cell(8,1);
timeCenterBurst_accumulator=cell(8,1);
ampBurst_accumulator=cell(8,1);
length_injected=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]; 
% Parallel pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'))
end
pool = parpool("Processes",8);
parfor jj=1:8
        [analogData,timeVals,analogData0,num_burst_gatherer,freq_burst_gatherer,...
        time_center_burst_gatherer,amp_burst_gatherer] = generate_synthetic_data_LFP(burstLen);
        [length_gatherer,length_gatherer_tf_l,length_gatherer_tf_h,...
         length_gatherer_MP,length_gatherer_MP_tf_l,length_gatherer_MP_tf_h,...
         length_gatherer_hilbert,length_gatherer_hilbert_tf_l,length_gatherer_hilbert_tf_h,...
         length_gatherer_wavelet,length_gatherer_wavelet_tf_l,length_gatherer_wavelet_tf_h] = run_all_methods(analogData,timeVals);
        %%%%%%%%%%%%%%%%%%%%%%%%%% OMP-GEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        length_accumulator_omp_gear{jj,1} =  length_gatherer; % default threshold used for analysis
        length_accumulator_omp_gear_tf_l{jj,1} = length_gatherer_tf_l; 
        length_accumulator_omp_gear_tf_h{jj,1} = length_gatherer_tf_h; 
        gabor_accmulator_omp_gear{jj,1}=gaborInfo_current;
        %%%%%%%%%%%%%%%% MP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        length_accumulator_MP{jj,1} = length_gatherer_MP; % default threshold used for analysis
        length_accumulator_MP_tf_l{jj,1} = length_gatherer_MP_tf_l;
        length_accumulator_MP_tf_h{jj,1} = length_gatherer_MP_tf_h;
        gabor_accmulator_MP{jj,1}=gaborInfo_MP;
        %%%%%%%%%%%%%%% Hilbert %%%%%%%%%%%%%%%%%%%%%%%%%%%
        length_accumulator_hilbert{jj,1} = length_gatherer_hilbert; % default threshold used for analysis
        length_accumulator_hilbert_tf_l{jj,1} = length_gatherer_hilbert_tf_l;
        length_accumulator_hilbert_tf_h{jj,1} = length_gatherer_hilbert_tf_h;
        %%%%%%%%%%%%%%% Wavelet %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % we do not use wavelet and feingold in our analysis, although the
        % result is consistent with these methods as well.
        length_accumulator_wavelet{jj,1} = length_gatherer_wavelet; % default threshold used for analysis
        length_accumulator_wavelet_tf_l{jj,1} = length_gatherer_wavelet_tf_l;
        length_accumulator_wavelet_tf_h{jj,1} = length_gatherer_wavelet_tf_h;
        %%%%%%%%%%%%%%%% Analong data- Synthetic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        analogData_accumulator{jj,1}=analogData;
        analogData0_accumulator{jj,1}=analogData0;
        timeVals_accumulator{jj,1}=timeVals;
        numBurst_accumulator{jj,1}=num_burst_gatherer;
        freqBurst_accumulator{jj,1}=freq_burst_gatherer;
        timeCenterBurst_accumulator{jj,1}=time_center_burst_gatherer;
        ampBurst_accumulator{jj,1}=amp_burst_gatherer;
end
save('LFP_synth_gamma_all_methods_elec41_all_ORI.mat','length_accumulator_omp_gear','length_accumulator_omp_gear_tf_l','length_accumulator_omp_gear_tf_h',...
    'gabor_accmulator_omp_gear','length_accumulator_MP','length_accumulator_MP_tf_l','length_accumulator_MP_tf_h',...
    'gabor_accmulator_MP','length_accumulator_hilbert','length_accumulator_hilbert_tf_l','length_accumulator_hilbert_tf_h',...
    'length_accumulator_wavelet','length_accumulator_wavelet_tf_l','length_accumulator_wavelet_tf_h',...
    'analogData_accumulator','analogData0_accumulator','timeVals_accumulator','numBurst_accumulator',...
    'freqBurst_accumulator','timeCenterBurst_accumulator','ampBurst_accumulator','length_injected','-v7.3');
