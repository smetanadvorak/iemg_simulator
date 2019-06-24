mf_mdl = MuscleForceMdl(N, rr, fs);

mf_mdl.plot_twitches();

fmax = mf_mdl.normalize_mvc(mvc_spikes);

f_offline = mf_mdl.generate_force_offline(lq_spikes);


%%
profile on

T = 25*fs;
excitation = linspace(0,1,T)';
spikes = zeros(T,N);
ipi = zeros(T,N);
prev_neuron_state = nan(N,1);

f_online = zeros(T,1);
mf_mdl.init_online_buffer();
for t = 1:T
    [spikes(t,:), prev_neuron_state, ipi(t,:)] = generate_spike_train_vect(t, prev_neuron_state, excitation(t), rt, minfr, maxfr, frs, CV, fs);
    f_online(t) = mf_mdl.generate_force_online(spikes(t,:), ipi(t,:));

    if ~mod(t,fs)
        fprintf('%d seconds generated\n', floor(t/fs));
    end
end

profile report