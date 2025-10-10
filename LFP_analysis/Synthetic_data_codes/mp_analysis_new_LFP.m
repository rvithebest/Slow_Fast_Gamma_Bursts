length_accumulator_omp_gear = cell(8,4);
gabor_accmulator_omp_gear=cell(8,1);
length_accumulator_MP=cell(8,4);
gabor_accmulator_MP=cell(8,1);
length_accumulator_hilbert=cell(8,4);
length_accumulator_wavelet=cell(8,4);
length_accumulator_feingold=cell(8,4);
analogData_accumulator=cell(8,3);
length_injected=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4];
for v = 1 % only once first 
    for jj=1:8
        clearvars -except v length_accumulator_omp_gear gabor_accmulator_omp_gear length_injected jj length_accumulator_MP gabor_accmulator_MP length_accumulator_hilbert analogData_accumulator length_accumulator_wavelet length_accumulator_feingold;
        % [analogData, timeVals, analogData0] = synthetic_data_LFP_new(length_injected(jj));
        % load the following file- "synth_data_length_pos_jj.mat"
        load(['synth_data_length_pos_', num2str(jj), '.mat']);
        MP_length_all_trials;   
        %%%%%%%%%%%%%%%%%%%%%%%%%% OMP-GEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        length_accumulator_omp_gear{jj,1} = [length_accumulator_omp_gear{jj,1}, length_gatherer];
        length_accumulator_omp_gear{jj,2} = [length_accumulator_omp_gear{jj,2}, length_gatherer_beta];
        length_accumulator_omp_gear{jj,3} = [length_accumulator_omp_gear{jj,3}, length_gatherer_gamma];
        length_accumulator_omp_gear{jj,4} = [length_accumulator_omp_gear{jj,4}, length_gatherer_delta];
        gabor_accmulator_omp_gear{jj,1}=gaborInfo_current;
        %%%%%%%%%%%%%%%% MP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        length_accumulator_MP{jj,1} = [length_accumulator_MP{jj,1}, length_gatherer_MP];
        length_accumulator_MP{jj,2} = [length_accumulator_MP{jj,2}, length_gatherer_MP_beta];
        length_accumulator_MP{jj,3} = [length_accumulator_MP{jj,3}, length_gatherer_MP_gamma];
        length_accumulator_MP{jj,4} = [length_accumulator_MP{jj,4}, length_gatherer_MP_delta];
        gabor_accmulator_MP{jj,1}=gaborInfo_MP;
        %%%%%%%%%%%%%%% Hilbert %%%%%%%%%%%%%%%%%%%%%%%%%%%
        length_accumulator_hilbert{jj,1} = [length_accumulator_hilbert{jj,1}, length_gatherer_hilbert];
        length_accumulator_hilbert{jj,2} = [length_accumulator_hilbert{jj,2}, length_gatherer_hilbert_beta];
        length_accumulator_hilbert{jj,3} = [length_accumulator_hilbert{jj,3}, length_gatherer_hilbert_gamma];
        length_accumulator_hilbert{jj,4} = [length_accumulator_hilbert{jj,4}, length_gatherer_hilbert_delta];
        %%%%%%%%%%%%%%% Wavelet %%%%%%%%%%%%%%%%%%%%%%%%%%%
        length_accumulator_wavelet{jj,1} = [length_accumulator_wavelet{jj,1}, length_gatherer_wavelet];
        length_accumulator_wavelet{jj,2} = [length_accumulator_wavelet{jj,2}, length_gatherer_wavelet_beta];
        length_accumulator_wavelet{jj,3} = [length_accumulator_wavelet{jj,3}, length_gatherer_wavelet_gamma];
        length_accumulator_wavelet{jj,4} = [length_accumulator_wavelet{jj,4}, length_gatherer_wavelet_delta];
        %%%%%%%%%%%%%%% Feingold %%%%%%%%%%%%%%%%%%%%%%%%%%%
        length_accumulator_feingold{jj,1} = [length_accumulator_feingold{jj,1}, length_gatherer_feingold];
        length_accumulator_feingold{jj,2} = [length_accumulator_feingold{jj,2}, length_gatherer_feingold_beta];
        length_accumulator_feingold{jj,3} = [length_accumulator_feingold{jj,3}, length_gatherer_feingold_gamma];
        length_accumulator_feingold{jj,4} = [length_accumulator_feingold{jj,4}, length_gatherer_feingold_delta];
        %%%%%%%%%%%%%%%% Analong data- Synthetic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        analogData_accumulator{jj,1}=[analogData_accumulator{jj,1},analogData];
        analogData_accumulator{jj,2}=[analogData_accumulator{jj,2},analogData0];
        analogData_accumulator{jj,3}=[analogData_accumulator{jj,3},timeVals];
        save('LFP_synth_gamma_all_methods_elec41_all_ORI.mat','length_accumulator_omp_gear','v','gabor_accmulator_omp_gear','length_injected','length_accumulator_MP','gabor_accmulator_MP','length_accumulator_hilbert','analogData_accumulator','jj','length_accumulator_wavelet','length_accumulator_feingold');
    end
end
