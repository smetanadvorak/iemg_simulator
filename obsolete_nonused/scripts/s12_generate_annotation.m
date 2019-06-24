%% Generate full causal annotation
% Problem here is that original spikes correspond to the neuron's spikes,
% while decomposed spikes usually correspond to some center of the action
% potential. Thus, decomposed spikes are around +5ms, compared to the
% neuron's spikes. This delay comes from neuron's axon and mf conduction
% velocities. In multichannel recordings, while axon delay is constant, mf
% conduction varies depending on the electrode's longitudinal position. I
% decided to have the same delay for all channels, averaging the mf
% conduction delay. 

% annotation_delays = zeros(N,1); % In samples!
% for m = 1:N
%     mf_conduction_delay = round(fs * abs(mean((MUs(m).nmj_z - mean(electrode_pts(:,3)))./MUs(m).mf_cv)));
%     axon_delay = round(fs * mean(MUs(m).mnap_delays));
%     annotation_delays(m) = mf_conduction_delay + axon_delay;
% end
% 
% annotation_spikes_full = nan(size(spikes));
% for m = 1:N
%     annotation_spikes_full(:,m) = shift_padding(spikes(:,m), annotation_delays(m), 1);
% end
% annotation_firings_full = spikes2firings(annotation_spikes_full);


%% Generate causal annotation for detectable MUs
% annotation_spikes_detectable = nan(size(spikes,1), numel(detectable_ind));
% for m = 1:numel(detectable_ind)
%     annotation_spikes_detectable(:,m) = shift_padding(spikes(:,detectable_ind(m)), annotation_delays(detectable_ind(m)), 1);
% end
% annotation_firings_detectable = spikes2firings(annotation_spikes_detectable);


%% Generate full ap-centered annotation
ap_centers_full = zeros(N, n_channels);
for m = 1:N
    temp_traj = traj_mixing_mat(0, traj_n_nodes, n_channels);
    ap_centers_full(m, :) = find_ap_center_gauss(MUs(m).muap * temp_traj' * electrode_diff_mat')';
end
spikes_centered_full = zeros([size(spikes), n_channels]);
for m = 1:N
    for ch = 1:n_channels
        spikes_centered_full(:,m,ch) = [zeros(ap_centers_full(m, ch)-1, 1); spikes(1:end-ap_centers_full(m, ch)+1, m)];
    end
end
firings_centered_full = spikes2firings(spikes_centered_full);

%% Generate ap-centered annotation for detectable motor units
spikes_centered_detectable = spikes_centered_full(:, detectable_ind, :);
firings_centered_detectable = spikes2firings(spikes_centered_detectable);


%% Plot EMG with spikes
figure; set(gcf, 'position', [70, 100, 1300, 600]);

separator_step = mean(std(emg)) * 5;
separator = (1:n_channels) * separator_step;
separator = repmat(separator, size(emg,1), 1);

% Plot EMG (separator for multichannel view)
plot(profile_timeline, emg + separator, 'linewidth', 1.1); hold on;
plot(profile_timeline, squeeze(sum(spikes_centered_detectable,2)) * std(emg(:)) + separator, 'linewidth', 1.2);
title(sprintf('Simulated EMG, [V], SNR=%2.2f', SNR));
xlabel('Time, s'); ylabel('Amplitude, arbitrary units');
legend('Signal', 'Centered Firings');

% APs labelling
% text_level = separator(1,end) + mean(std(emg)) * 5;
% for t = 1:size(spikes,1)
%     for m = 1:numel(detectable_ind)
%         if annotation_spikes_detectable(t,m)
%             text(profile_timeline(t), text_level + m*mean(std(emg))/4, num2str(m));
%         end
%     end
% end
clear separator t m;

%% Plot dictionary across trajectory



