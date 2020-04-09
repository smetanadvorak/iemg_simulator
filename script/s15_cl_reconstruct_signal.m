temp_reconstructed_emg_full = zeros(size(spikes,1), electrode.n_points * electrode.n_nodes);
temp_reconstructed_emg_detectable = zeros(size(spikes,1), electrode.n_points * electrode.n_nodes);

for mu = 1:mu_pool.N
    for pt = 1:electrode.n_points * electrode.n_nodes
        temp_train = conv(spikes(:,mu), MUs(mu).muap(:,pt), 'full');
        temp_reconstructed_emg_full(:,pt) = temp_reconstructed_emg_full(:,pt) + temp_train(1: end-size(MUs(mu).muap, 1)+1);                  
        if any(mu == detectable_ind)
            temp_reconstructed_emg_detectable(:,pt) = temp_reconstructed_emg_detectable(:,pt) + temp_train(1: end-size(MUs(mu).muap, 1)+1);                  
        end
    end
end

% reconstructed_emg_full = reconstructed_emg_full * electrode.diff_mat';
% reconstructed_emg = reconstructed_emg * electrode.diff_mat';
reconstructed_emg_full = zeros(size(temp_reconstructed_emg_full,1), electrode.n_channels);
reconstructed_emg_detectable = zeros(size(temp_reconstructed_emg_detectable,1), electrode.n_channels);

traj_parameter_map = force;
for i = 1:T
    temp_traj = electrode.traj_mixing_mat(traj_parameter_map(i), electrode.n_nodes, electrode.n_channels);
    reconstructed_emg_full(i,:) = temp_reconstructed_emg_full(i,:) * temp_traj' * electrode.diff_mat';
    reconstructed_emg_detectable(i,:) = temp_reconstructed_emg_detectable(i,:) * temp_traj' * electrode.diff_mat';
end



% figure; set(gcf, 'position', [70, 100, 1300, 600]);
% 
% separator_step = mean(std(reconstructed_emg_full)) * 5;
% separator = (1:electrode.n_channels) * separator_step;
% separator = repmat(separator, size(reconstructed_emg_full,1), 1);
% 
% % Plot EMG (separator for multichannel view)
% subplot(2,1,1);
% plot(profile.timeline, emg + separator, 'k', 'linewidth',1); hold on;
% plot(profile.timeline, reconstructed_emg_detectable + separator, 'linewidth', 1, 'color', 'b'); hold on;
% legend('Original signal', 'Reconstr. from detectable dict.');
% 
% subplot(2,1,2);
% plot(profile.timeline, emg + separator, 'k', 'linewidth',1); hold on;
% plot(profile.timeline, reconstructed_emg_full + separator, 'linewidth', 1, 'color', 'g'); hold on;
% legend('Original signal', 'Reconstr. from full dict.');
% 
% title(sprintf('Reconstructed EMG vs Simulated'));
% xlabel('Time, s'); ylabel('Amplitude, Arbitrary units');



clear temp* reconstructed_* separator* mu
