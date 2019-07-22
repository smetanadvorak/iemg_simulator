load('/Users/akmbpro/Nextcloud/EMG Data/simulated/physiological/Paper_On_Simulation/single_channel/single_channel_mdl.mat');

targetMU = 28;
to_plot = 1/1000*fs : 7/1000*fs;
timeline = (0:numel(MUs(targetMU).calc_muap(0))-1)/fs*1000;
muap_original = MUs(targetMU).calc_muap(0) * electrode.diff_mat';

figure; hold on;
% Plot jittered versions
for i = 1:10
    muap_to_plot = MUs(targetMU).calc_muap(35e-6) * electrode.diff_mat';
    p2=plot(timeline(to_plot), muap_to_plot(to_plot, :), 'b', 'linewidth', 0.2);
end

% Plot original MUAP
p1=plot(timeline(to_plot), muap_original(to_plot, :) , 'k', 'linewidth', 1);

legend([p1,p2],'Original MUAP',  '10 APs with jitter');
xlabel('Time, ms'); ylabel('Amplitude, arbitrary units');
axis tight


set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', [0, 0, 10.5, 7.4]); %[0, 0, 7.4, 5.2], [0, 0, 10.5, 7.4]);
set(gcf, 'papersize', [21 29.7]);
set(gca, 'fontunits', 'points', 'fontsize', 8);
print -depsc2 jitter_std.eps


%%
figure; hold on;
% Create jittered versions
n_jitters = 100;
jittered_muaps = zeros(size(muap_original,1), n_jitters);
for i = 1:n_jitters
    jittered_muaps(:,i) = MUs(targetMU).calc_muap(35e-6) * electrode.diff_mat';
end

% Plot std
p2 = plot(timeline(to_plot), muap_original(to_plot) + 3*std(jittered_muaps(to_plot,:), [], 2), 'b');
plot(timeline(to_plot), muap_original(to_plot) - 3*std(jittered_muaps(to_plot,:), [], 2), 'b','HandleVisibility','off');

% Plot original MUAP
p1 = plot(timeline(to_plot), muap_original(to_plot, :) , 'k', 'linewidth', 1);

legend([p1, p2], 'Original MUAP', '$\pm 3 \sigma$ interval');
xlabel('Time, ms'); ylabel('Amplitude, arbitrary units');
axis tight


set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', [0, 0, 10.5, 7.4]); %[0, 0, 7.4, 5.2], [0, 0, 10.5, 7.4]);
set(gcf, 'papersize', [21 29.7]);
set(gca, 'fontunits', 'points', 'fontsize', 8);
print -depsc2 jitterred_muaps.eps

clear  muap_to_plot timeline n_to_plot