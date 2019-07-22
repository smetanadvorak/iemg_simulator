%load('/Users/akmbpro/Nextcloud/EMG Data/simulated/physiological/Paper_On_Simulation/territories_assessement/multichannel_scan_mdl.mat')

%% Upsample trajectory nodes
upsample = 1;
snapshots = 0 : 1/numel(electrode.traj_transforms(:,1))/upsample : 1;
spatial_step = electrode.traj_step/upsample;
n_triggers = 30;

%% Estimate diameters
mus_to_take = round(2*mu_pool.N/3);
estimated_diameters = nan(mu_pool.N,1);

for m = 1:mus_to_take
    % Calculate MUAP in all positions (including upsampled positions)
    muap_z = zeros(size(MUs(m).muap, 1), length(snapshots));
    for i = 1 : length(snapshots)
        muap_z(:,i) = MUs(m).muap * electrode.traj_mixing_mat(snapshots(i), electrode.n_nodes, electrode.n_channels)' * electrode.diff_mat';
    end
    
    % Calculate std of the MUAP in each position
    %muap_std = std(muap_z);
    
    muap_std = zeros(length(snapshots),1);
    for n = 1:length(snapshots)
        valid_part = muap_z(abs(muap_z(:,n)) > emg_noise_std/sqrt(n_triggers),n);
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


%% Linear fit
color_scale_min = min(abs(mu_pool.mn_pool.centers(:,2))/Rmuscle);
color_scale_max = max(abs(mu_pool.mn_pool.centers(:,2))/Rmuscle - color_scale_min);
k = true_diameters(~isnan(estimated_diameters))\estimated_diameters(~isnan(estimated_diameters));
figure; hold on
plot(true_diameters, estimated_diameters, '.');
for i = 1:mus_to_take
    col = ((abs(mu_pool.mn_pool.centers(i,2)))/Rmuscle - color_scale_min)/color_scale_max;
    text(true_diameters(i)-0.1, estimated_diameters(i), num2str(i), 'fontsize', 12, 'color', [0, col, 1-col]);
   % plot(true_diameters(i), estimated_diameters(i), 'o', 'markersize', 20, 'color', [0, col, 1-col]);
end
hp = plot([0;true_diameters], k*[0;true_diameters], 'k');
title(''); 
xlabel('Diameters estimated from exact MF coordinates, mm'); 
ylabel('Diameters estimated using scanning electrode, mm');
legend(hp, sprintf('Linear fit: k=%2.2f', k));


%%
figure;
for i = 1:mus_to_take
    plot(i, estimated_diameters(i), 'ko'); hold on;
    plot(i, true_diameters(i), 'kx');
    line([i, i], [estimated_diameters(i), true_diameters(i)], 'linestyle','--');
end
legend('Diameters estimated from scanning', 'Model diameters $\left(2\sqrt{a_i/\pi}\right)$');
title('Diameters of detectable motor units');
xlabel('Motor unit');
ylabel('Diameter, mm');

%% Figure with crossed territories;


%%
mu_pool.show_innervation_areas_2d(find(~isnan(estimated_diameters)));
plot(electrode.pts(1:end, 1), electrode.pts(:,2), '*k', 'linewidth', 1);
%%
%clear upsample snapshots m muap_z first_detection last_detection

