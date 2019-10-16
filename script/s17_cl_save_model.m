
cd 'simulation_output';

% save('simulation_mdl', 'fs', 'fsl', 'dt', 'dz');
% save('simulation_mdl', 'mu_pool', '-append');
% 
% save('simulation_mdl', 'MUs', 'max_muap_len', '-append');
% save('simulation_mdl', 'Rmuscle', 'Lmuscle', '-append');
% 
% save('simulation_mdl', 'mf_mdl', '-append');
% save('simulation_mdl', 'e2fModel', 'pidc', '-append');
% save('simulation_mdl', 'electrode', '-append');
%                    
% save('simulation_mdl', 'mvc_emg_std', '-append');

save('simulation_mdl', 'fs', 'fsl', 'dt', 'dz', ...
                       'mu_pool', 'MUs', 'max_muap_len', ...
                       'Rmuscle', 'Lmuscle', ...
                       'mf_mdl', 'e2fModel', 'pidc', ...
                       'electrode','mvc_emg_std','emg_noise_std', 'SNR');

cd ..