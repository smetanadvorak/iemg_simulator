%% Generate full centered and shortened dictionary
if n_channels == 1
    
    dict_muap_half_len = floor(10/2/1000 * fs);
    dict_muap_full_len = 2*dict_muap_half_len + 1;
    
    dictionary_for_decomp_full = zeros(dict_muap_full_len, N);
    dictionary_for_decomp_detectable = zeros(dict_muap_full_len, numel(detectable_ind));
    
    ap_centers_full = zeros(N,1);
    for m = 1:N
        temp_traj = traj_mixing_mat(0, traj_n_nodes, n_channels);
        ap_centers_full(m) = find_ap_center_gauss(dictionary_full_init(:,m));
        dictionary_for_decomp_full(:,m) = take_scope(dictionary_full_init(:,m), ap_centers_full(m), dict_muap_half_len);
    end
    
    ap_centers_detectable = zeros(numel(detectable_ind),1);
    for m = 1:size(detectable_ind,1)
        temp_traj = traj_mixing_mat(0, traj_n_nodes, n_channels);
        ap_centers_detectable(m) = find_ap_center_gauss(dictionary_detectable_init(:,m));
        dictionary_for_decomp_detectable(:,m) = take_scope(dictionary_detectable_init(:,m), ap_centers_detectable(m), dict_muap_half_len);
    end
    
    % Generate annotation
    firings_for_decomp_full = shift_firings(spikes2firings(spikes), - dict_muap_half_len + ap_centers_full);
    firings_for_decomp_detectable = shift_firings(spikes2firings(spikes(:,detectable_ind)), - dict_muap_half_len + ap_centers_detectable);
    spikes_for_decomp_full = firings2spikes(firings_for_decomp_full, length(spikes));
    spikes_for_decomp_detectable = firings2spikes(firings_for_decomp_detectable, length(spikes));
    
    
    %% Reconstruct it
    reconstructed_for_decomp_full = conv_trains_muaps(spikes_for_decomp_full, dictionary_for_decomp_full, 'causal');
    reconstructed_for_decomp_detectable = conv_trains_muaps(spikes_for_decomp_detectable, dictionary_for_decomp_detectable, 'causal');
    
     
    
    figure; set(gcf, 'position', [70, 100, 1300, 600]);
    
    separator_step = mean(std(reconstructed_for_decomp_full)) * 5;
    separator = (1:n_channels) * separator_step;
    separator = repmat(separator, size(reconstructed_for_decomp_full,1), 1);
    
    % Plot EMG (separator for multichannel view)
    subplot(2,1,1);
    plot(profile_timeline, emg + separator, 'k', 'linewidth',1); hold on;
    plot(profile_timeline, reconstructed_for_decomp_detectable + separator, 'linewidth', 1, 'color', 'b'); hold on;
    plot(profile_timeline, squeeze(sum(spikes_for_decomp_detectable,2)) * std(emg(:)) + separator, 'linewidth', 1.2);
    legend('Original signal', 'Reconstr. from detectable dict.');
    
    subplot(2,1,2);
    plot(profile_timeline, emg + separator, 'k', 'linewidth',1); hold on;
    plot(profile_timeline, reconstructed_for_decomp_full + separator, 'linewidth', 1, 'color', 'g'); hold on;
    plot(profile_timeline, squeeze(sum(spikes_for_decomp_full,2)) * std(emg(:)) + separator, 'linewidth', 1.2);
    legend('Original signal', 'Reconstr. from full dict.');
    
    title(sprintf('Reconstructed EMG vs Simulated'));
    xlabel('Time, s'); ylabel('Amplitude, Arbitrary units');
    
    clear temp*
    
    
end
