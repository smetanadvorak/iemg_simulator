%% Select detectable MUAPs
% Annotation should contain only those MUs, who's APs are prominent in the
% signal and, thus, decomposable. In fact, prominence is a tricky notion,
% since small MU's AP are prominent at their recruitment level and are not
% prominent at higher excitation due to appearance of larger MUs;
% So, prominence of a MUAP should be:
% 1) Assessed at its MU recruitment threshold;
% 2) Compared to the instrumentation noise. 

% First criterion: by explained variance of the signal at the recruitment
% threshold

% Power of motor unit's signal is calculated as it's
% MUAP's energy times the firing rate at recruitment.
% Contribution is then calculated as percentage of that power in the total
% power of the emg signal at that excitation level.
total_explained_variance = 0.99;
over_noise = 6;

recr_contribution = zeros(N,1);
expl_detectable = zeros(N,1);
expl_detectable(1) = 1;

for m = 2:N
    fr = zeros(m,1);
	explained_variance = zeros(m,n_channels);
    % Get firing rates and powers of all MUs that are active at this 
    % excitation level (basically, all smaller MUs)
    for i = 1:m
        fr(i) = calculate_fr(rt(i), minfr(i), maxfr(i), frs(i), rt(m) + eps); % excitation should be > rt to activate the motor unit
        muap_sum_squares = sum((MUs(i).muap * eye(size(traj_mixing_mat(0, traj_n_nodes, n_channels)')) * electrode_diff_mat').^2);
        explained_variance(i,:) = muap_sum_squares .* fr(i);
    end
    explained_variance = explained_variance./repmat(sum(explained_variance), m, 1); % Normalize contributions 
%    [explained_variance, ind] = sort(explained_variance, 'descend');
%    explained_variance = cumsum(explained_variance);
%    recr_detectable(m) = any(explained_variance(ind(m),:) < total_explained_variance);
    expl_detectable(m) = any(explained_variance(m,:) > 1-total_explained_variance);
end

% By prominence compared to instrumentation noise:
prom_detectable = zeros(N,1);
for m = 1:N
    prom_detectable(m) = any( max(abs((MUs(m).muap * electrode_diff_mat'))) > over_noise * emg_noise_std );
    %prom_detectable(m) = any( sqrt(mean( (MUs(m).muap * electrode_diff_mat').^2 )) > emg_noise_std );
end

detectable = prom_detectable .* expl_detectable;
detectable_ind = find(detectable);


