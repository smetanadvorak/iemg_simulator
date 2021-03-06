
%% Upsample trajectory nodes
upsample = 5;
snapshots = 0 : 1/numel(electrode.traj_transforms(:,1))/upsample : 1;
spatial_step = electrode.traj_step/upsample;

%% Estimate diameters
estimated_diameters = zeros(numel(prom_detectable_ind),1);

for m = 1:numel(prom_detectable_ind)
    % Calculate MUAP in all positions (including upsampled positions)
    muap_z = zeros(size(MUs(prom_detectable_ind(m)).muap, 1), length(snapshots));
    for i = 1 : length(snapshots)
        muap_z(:,i) = MUs(prom_detectable_ind(m)).muap * electrode.traj_mixing_mat(snapshots(i), electrode.n_nodes, electrode.n_channels)' * electrode.diff_mat';
    end
    
    % Calculate std of the MUAP in each position
    %muap_std = std(muap_z);
    
    muap_std = zeros(length(snapshots),1);
    for n = 1:length(snapshots)
        valid_part = muap_z(abs(muap_z(:,n)) > emg_noise_std,n);
        if isempty(valid_part) 
            muap_std(n) = 0;
        else
            muap_std(n) = std( valid_part );
        end
    end
    
    % Find the first and last positions where MUAP is detectable and above threshold
    %first_detection = find( (muap_std > max(muap_std)/50) & (muap_std > 0*emg_noise_std), 1, 'first');
    %last_detection = find( (muap_std > max(muap_std)/50) & (muap_std > 0*emg_noise_std), 1, 'last');

    first_detection = find( max(abs(muap_z))' > 4*emg_noise_std & (muap_std > 1*emg_noise_std), 1, 'first');
    last_detection = find(  max(abs(muap_z))' > 4*emg_noise_std & (muap_std > 1*emg_noise_std), 1, 'last');
    
    if isempty(first_detection) || isempty(last_detection)
        estimated_diameters(m) = nan;
    else
        estimated_diameters(m) = (last_detection - first_detection) * spatial_step;
    end
end

%% Plot estimated diameters and diameters calculated from true areas
figure; 
%plot(prom_detectable_ind, estimated_diameter, 'o', 'markersize', 5); hold on; 
plot(estimated_diameters, 'linewidth', 2); hold on; 
true_areas = 2*sqrt(mu_pool.calc_innervation_areas_res('polygone_area'));
%true_areas = 2*sqrt(mu_pool.innervation_areas);
true_diameters = true_areas(prom_detectable_ind);

%plot(prom_detectable_ind, res_areas(prom_detectable_ind), 'o', 'markersize', 5);
plot(true_diameters, 'linewidth', 2);
legend('Diameters estimated from scanning', 'Model diameters $\left(2\sqrt{a_i/\pi}\right)$');
title('Diameters of detectable motor units');

%% Linear fit
color_scale_min = min(abs(mu_pool.mn_pool.centers(:,2))/Rmuscle);
color_scale_max = max(abs(mu_pool.mn_pool.centers(:,2))/Rmuscle - color_scale_min);
k = true_diameters(~isnan(estimated_diameters))\estimated_diameters(~isnan(estimated_diameters));
figure; hold on
plot(true_diameters, estimated_diameters, '.');
for i = 1:numel(prom_detectable_ind)
    col = ((abs(mu_pool.mn_pool.centers(prom_detectable_ind(i),2)))/Rmuscle - color_scale_min)/color_scale_max;
    text(true_diameters(i)-0.1, estimated_diameters(i), num2str(prom_detectable_ind(i)), 'fontsize', 12, 'color', [0, col, 1-col]);
   % plot(true_diameters(i), estimated_diameters(i), 'o', 'markersize', 20, 'color', [0, col, 1-col]);
end
hp = plot([0;true_diameters], k*[0;true_diameters], 'k');
title(''); 
xlabel('Diameters estimated from exact MF coordinates, mm'); 
ylabel('Diameters estimated using scanning electrode, mm');
legend(hp, sprintf('Linear fit: k=%2.2f', k));


%%
figure;
for i = 1:numel(prom_detectable_ind)
    plot(prom_detectable_ind(i), estimated_diameters(i), 'ko'); hold on;
    plot(prom_detectable_ind(i), true_diameters(i), 'kx');
    line([prom_detectable_ind(i), prom_detectable_ind(i)], [estimated_diameters(i), true_diameters(i)], 'linestyle','--');
end
legend('Diameters estimated from scanning', 'Model diameters $\left(2\sqrt{a_i/\pi}\right)$');
title('Diameters of detectable motor units');
xlabel('Motor unit');
ylabel('Diameter, mm');

%% Figure with crossed territories;


%%
mu_pool.show_innervation_areas_2d(prom_detectable_ind);
%%
clear upsample snapshots m muap_z first_detection last_detection

