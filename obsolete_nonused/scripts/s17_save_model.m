
cd 'simulation_output';

save('simulation_mdl', 'fs', 'dt', 'dz', 'CV');
save('simulation_mdl', 'N', 'rr', 'rm', 'sz', 'rt', 'minfr', 'maxfr', 'frs', 'innerv_area', '-append');

save('simulation_mdl', 'Rmuscle', 'Lmuscle', 'Nmf', 'innerv_num', '-append');
save('simulation_mdl', 'MUs', 'max_muap_len', '-append');
save('simulation_mdl', 'Tf', 'Pf', 'twitches_mat', 'twitches_cell', 'fmax', '-append');
save('simulation_mdl', 'e2fModel', 'pidc', '-append');
save('simulation_mdl', 'electrode_type', 'electrode_step', 'electrode_pts_init', 'simplified_pts', ...
                       'electrode_diff_mat', 'n_electrodes', 'electrode_diff_mat', 'n_channels', ...
                       'n_points', 'traj_n_nodes', 'traj_transforms', 'traj_mixing_fun',...
                       'traj_mixing_mat', 'electrode_type_short', '-append');
                   
save('simulation_mdl', 'mvc_emg_std', '-append');
%save('simulation_mdl', 'T', 'profile_timeline', 'force_profile', 'profile_type', 'traj_parameter_map', ...
%                       'traj_parameter_map_type', '-append');
                   
save('simulation_mdl', 'detectable_ind', '-append');                   




cd ..