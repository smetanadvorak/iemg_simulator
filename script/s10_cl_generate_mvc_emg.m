%% Generate MVC spikes
mvc_spikes = mu_pool.mn_pool.generate_spike_train_gauss(1:(5*fs), nan(N,1), ones(5*fs,1), fs);

%% Generate MVC EMG
mvc_emg = zeros(size(mvc_spikes,1), electrode.n_points * electrode.n_nodes);
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
mvc_emg = mvc_emg * electrode.traj_mixing_mat(0.5, electrode.n_nodes, electrode.n_channels)' * electrode.diff_mat';
mvc_emg_std = std(mvc_emg(fs:end-fs, :));

%% Noise
SNR = 13;
emg_noise_std = mvc_emg_std * 10^(-SNR/10);
emg_noise_std = mean(emg_noise_std) * ones(size(emg_noise_std)); % Average in all channels


clear i t muap_to_add