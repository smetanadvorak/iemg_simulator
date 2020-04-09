nmj_jitter = 0; %35e-6;

%% Generate EMG according to profile
emg = zeros(size(spikes,1), electrode.n_points * electrode.n_nodes);
T = size(spikes,1);
traj_parameter_map = linspace(0,1,T);

for i = 1:mu_pool.N
    for t = 1:T - max_muap_len - 1
        if spikes(t,i)
            muap_to_add = MUs(i).calc_muap(nmj_jitter);
            emg(t:t+size(muap_to_add,1)-1, :) = ...
                    emg(t:t+size(muap_to_add,1)-1, :) + muap_to_add;
        end
    end
    fprintf('%d MUAP trains are generated\n', i);
end

emg_clear = zeros(size(emg,1), electrode.n_channels);
for i = 1:T
    emg_clear(i,:) = emg(i,:) * electrode.traj_mixing_mat(traj_parameter_map(i), electrode.n_nodes, electrode.n_channels)' * electrode.diff_mat';
end

emg_noise = repmat(emg_noise_std, T, 1) .* randn(size(emg_clear));
emg = emg_clear + emg_noise;

clear muap_to_add i t emg_noise




