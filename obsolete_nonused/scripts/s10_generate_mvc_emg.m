

%% Generate MVC EMG
mvc_emg = zeros(size(mvc_spikes,1), n_points * traj_n_nodes);
for i = 1:N
    for t = 1:size(mvc_spikes,1) - max_muap_len - 1
        if mvc_spikes(t,i)
            muap_to_add = MUs(i).calc_muap(0);
            mvc_emg(t:t+size(muap_to_add,1)-1, :) = ...
                    mvc_emg(t:t+size(muap_to_add,1)-1, :) + muap_to_add;
        end
    end
    fprintf('%d MUAP trains are generated\n', i);
end

%%
mvc_emg = mvc_emg * traj_mixing_mat(0.5, traj_n_nodes, n_channels)' * electrode_diff_mat';
mvc_emg_std = std(mvc_emg(fs:end-fs, :));

clear i t muap_to_add