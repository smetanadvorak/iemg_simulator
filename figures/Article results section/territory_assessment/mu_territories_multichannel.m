%load('/Users/akmbpro/Nextcloud/EMG Data/simulated/physiological/Paper_On_Simulation/territories_assessement/multichannel_1mm.mat')

%% Upsample trajectory nodes
upsample = 1;
spatial_step = 1;
n_triggers = 1;
SNR = 10;
emg_noise_std = mvc_emg_std * 10^(-SNR/10);
emg_noise_std = mean(emg_noise_std);
mus_to_take = 1:round(2*mu_pool.N/3);
mus_to_exclude = [10, 23, 46, 54, 63, 67];
%electrode.diff_mat = [zeros(electrode.n_channels,1), eye(electrode.n_channels)];
%% Estimate diameters

estimated_diameters = nan(mu_pool.N,1);
for m = mus_to_take
    % Calculate MUAP in all positions (including upsampled positions)
    muap_z = MUs(m).muap * electrode.traj_mixing_mat(0, electrode.n_nodes, electrode.n_channels)' * electrode.diff_mat';
   
    % Calculate std of the MUAP in each position
    muap_std = zeros(size(muap_z,2),1);
    for n = 1:size(muap_z,2)
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
    last_detection  =  find( max(abs(muap_z))' > 4*emg_noise_std & (muap_std > 1*emg_noise_std), 1, 'last');
    
    if isempty(first_detection) || isempty(last_detection)
        estimated_diameters(m) = nan;
    else
        estimated_diameters(m) = (last_detection - first_detection) * spatial_step;
    end
end


plot(estimated_diameters, 'linewidth', 2); hold on; 
true_areas = mu_pool.calc_innervation_areas_res('confidence_ellipse', 0.99); 
true_diameters = 2*sqrt(true_areas/pi);

estimated_diameters(mus_to_exclude) = nan;
estimated_diameters(estimated_diameters < 1) = 1;

%% Figure with crossed territories;
hps = mu_pool.show_innervation_areas_2d(find(~isnan(estimated_diameters)));
plot(electrode.pts(1:end, 1), electrode.pts(:,2), '*k', 'markersize', 2, 'linewidth', 0.25);
xlim([-5.5;5.5]); ylim([-3;3]);

set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', [0, 0, 10.5, 7.4]); %[0, 0, 7.4, 5.2], [0, 0, 10.5, 7.4]);
set(gcf, 'papersize', [21 29.7]);
set(gca, 'fontunits', 'points', 'fontsize', 8);
print -depsc2 territories_assessment_circle.eps


%%
figure;
k = 0;
for i = mus_to_take
    if ~isnan(estimated_diameters(i)) 
        k = k + 1;
        po = plot(i, estimated_diameters(i), 'o', 'markersize', 2, 'color', hps(k).Color); hold on;
        px = plot(i, true_diameters(i), 'kx', 'markersize', 4, 'color', hps(k).Color); 
        if estimated_diameters(i) > true_diameters(i)
            text(i, true_diameters(i)-0.125, num2str(i), 'fontsize', 6, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'color', hps(k).Color);
        else
            text(i, true_diameters(i)+0.125, num2str(i), 'fontsize', 6, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'color', hps(k).Color);
        end
        line([i, i], [estimated_diameters(i), true_diameters(i)], 'color', 'k', 'linestyle',':');
    end
end
legend([po, px], '$\hat D$', '$D=\left(2\sqrt{a_i/\pi}\right)$', 'location', 'nw');
%title('Diameters of detectable motor units');
xlabel('Motor unit');
ylabel('Diameter, mm');
grid on;

set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', [0, 0, 10.5, 7.4]); %[0, 0, 7.4, 5.2], [0, 0, 10.5, 7.4]);
set(gcf, 'papersize', [21 29.7]);
set(gca, 'fontunits', 'points', 'fontsize', 8);
print -depsc2 territories_assessment_bar.eps



%%
%clear upsample snapshots m muap_z first_detection last_detection

