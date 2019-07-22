figure; set(gcf, 'position', [70, 100, 1300, 600]);

t_to_show = 1*fs : 5*fs;
channels_to_show = [1:1:1];%2:1:5;


separator_step = mean(std(emg(t_to_show,channels_to_show))) * 12;
separator = (1:numel(channels_to_show)) * separator_step;
separator = repmat(separator, size(emg(t_to_show,channels_to_show),1), 1);

% Plot EMG (separator for multichannel view)
plot(profile_timeline(t_to_show) - t_to_show(1)/fs, emg(t_to_show,channels_to_show) + separator, 'linewidth', 1, 'linestyle', '-'); hold on;
%plot(profile_timeline(t_to_show) - t_to_show(1)/fs, squeeze(sum(spikes_centered_detectable(t_to_show,:,channels_to_show),2)) * mean(std(emg(t_to_show,channels_to_show))) + separator, 'k', 'linewidth', 0.75);
%title(sprintf('Simulated EMG, [V], SNR=%2.2f', SNR));
xlabel('Time, s'); ylabel('Amplitude, Arbitrary units');
%legend('Signal', 'Centered Firings');

    % set(gca, 'YTick', (1:numel(channels_to_show)) * separator_step);
% yticklabels = cell(numel(channels_to_show),1);
% for i = 1:numel(channels_to_show)
%     yticklabels{i} = ['Ch. ', num2str(i)];
% end
% set(gca, 'YTickLabel', yticklabels);

% APs labelling
% text_level = separator(1,end) + mean(std(emg)) * 5;
% for t = 1:size(spikes,1)
%     for m = 1:numel(detectable_ind)
%         if annotation_spikes_detectable(t,m)
%             text(profile_timeline(t), text_level + m*mean(std(emg))/4, num2str(m));
%         end
%     end
% end
axis tight
grid minor
clear separator t m;