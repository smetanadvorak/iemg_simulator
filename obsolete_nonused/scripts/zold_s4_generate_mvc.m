T = 10*fs;
mvc_excitation = ones(T,1);

%% Calculate spike trains
mvc_spikes = zeros(T,N);
for m = 1:N
    mvc_prev_state = nan;
    for t = 1:T
        [mvc_spikes(t,m), mvc_prev_state] = generate_spike_train( t, mvc_prev_state, mvc_excitation(t), rt(m), minfr(m), maxfr(m), frs(m), CV, fs );
    end
    fprintf('%d Step spike trains generated\n', m);
end
clear mvc_prev_state

%% IPI signal generation out of spikes signal (for gain nonlinearity)
[~, mvc_ipi] = sawtooth2ipi(spikes2sawtooth([mvc_spikes(2:end, :); zeros(1,N)]));
mvc_gain = ones(size(mvc_spikes));
for i = 1:N
    mvc_gain(:,i) = get_force_gain_fuglevand(mvc_ipi(:,i), Tf(i));
end
clear mvc_ipi

%% Generate force
mvc_force = zeros(T, 1);
for i = 1:N
    for t = 1:T% - max_twitch_len - 1
        if mvc_spikes(t,i)
            twitch_to_add = twitches_cell{i};
            to_take = min(length(twitch_to_add), T - t);
            mvc_force(t:(t+to_take-1))  = ...
                    mvc_force(t:(t+to_take-1)) + mvc_gain(t,i) * twitch_to_add(1:to_take);
        end
    end
    fprintf('%d Step twitch trains are generated\n', i);
end

clear to_take twitch_to_add i t

figure; 
subplot(2,1,1); plot(mvc_force); hold on; plot(mvc_excitation * mean(mvc_force(mvc_excitation>0))); title('Identification data');


%% Estimate the MVC for normalization:
fmax = mean(mvc_force(fs:end-fs));

