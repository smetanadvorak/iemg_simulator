
%% Create a muscle model object
mf_mdl = MuscleForceMdl(mu_pool.mn_pool.N, mu_pool.mn_pool.rr, fs);
mf_mdl.plot_twitches();


%% Estimate the maxumum voluntary contraction, normalize force
% Define MVC spikes
mvc_T = 10*fs;
mvc_excitation = ones(mvc_T,1);
profile on
mvc_spikes = mu_pool.mn_pool.generate_spike_train_gauss(1:mvc_T, nan(mu_pool.mn_pool.N,1), mvc_excitation, fs);
profile report
[mvc_force, ~] = mf_mdl.normalize_mvc(mvc_spikes);

figure; 
plot((1:mvc_T)/fs, mvc_force); hold on; title('MVC after normalization');
ylabel('Force, normalized'); xlabel('Time');


clear mvc_*
%% Get the quasistatic inverse model of the muscle
qsi_T = 25*fs;
qsi_excitation = linspace(0,1,qsi_T)';
qsi_spikes = mu_pool.mn_pool.generate_spike_train_gauss(1:qsi_T, nan(mu_pool.mn_pool.N,1), qsi_excitation, fs);
qsi_force = mf_mdl.generate_force_offline(qsi_spikes);

% Invert the curve, get polynomial interpolation that passes through zero
%qsi_model = polyfit(qsi_force, qsi_excitation, 3);
% qsi_model = [qsi_force.^4, qsi_force.^3, qsi_force.^2, qsi_force.^1]\qsi_excitation;
% qsi_model = [qsi_model; 0];

mf_mdl.init_quasistatic_e2f_f2e_models(mu_pool);

figure; %plot(qsi_force, qsi_excitation); hold on; 
%plot(qsi_force, polyval(qsi_model, qsi_force), 'linewidth', 2); hold on;
plot(qsi_force, mf_mdl.f2e(qsi_force), 'linewidth', 2); hold on;
plot(qsi_force, qsi_excitation, 'linewidth', 2);
axis tight
ylabel('Excitation');
xlabel('Force, normalized');
title('Force response to quasistatic linearly increasing excitation and its polynomial fit');

%clear qsi_force qsi_excitation qsi_T qsi_spikes