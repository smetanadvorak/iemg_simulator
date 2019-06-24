
%% Create a muscle model object
mf_mdl = MuscleForceMdl(N, rr, fs);
mf_mdl.plot_twitches();



%% Estimate the maxumum voluntary contraction, normalize force
% Define MVC spikes
mvc_T = 10*fs;
mvc_excitation = ones(mvc_T,1);
mvc_spikes = generate_spike_train_vect(1:mvc_T, nan(N,1), mvc_excitation, rt, minfr, maxfr, frs, CV, fs);
[mvc_force, fmax] = mf_mdl.normalize_mvc(mvc_spikes);

figure; 
plot((1:mvc_T)/fs, mvc_force); hold on; title('MVC after normalization');
ylabel('Force, normalized'); xlabel('Time');


clear mvc_* 
%% Get the quasistatic inverse model of the muscle
qsi_T = 25*fs;
qsi_excitation = linspace(0,1,qsi_T)';
qsi_spikes = generate_spike_train_vect(1:qsi_T, nan(N,1), qsi_excitation, rt, minfr, maxfr, frs, CV, fs);
qsi_force = mf_mdl.generate_force_offline(qsi_spikes);

% Invert the curve, get polynomial interpolation
qsi_model = polyfit(qsi_force, qsi_excitation, 3);

figure; plot(qsi_force, qsi_excitation); hold on; plot(qsi_force, polyval(qsi_model, qsi_force), 'linewidth', 2);
axis tight
ylabel('Excitation');
xlabel('Force, normalized');
title('Force response to quasistatic linearly increasing excitation and its polynomial fit');

clear qsi_force qsi_excitation qsi_T qsi_spikes