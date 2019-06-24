%% Define a slowly growing linear excitation
T = 25 * fs;
lq_exc = linspace(0,1,T)';


%% Calculate the spikes
lq_spikes = zeros(T,N); % Identification data
for m = 1:N
    prev_state = nan;
    for t = 1:T
        [lq_spikes(t,m), prev_state] = generate_spike_train( t, prev_state, lq_exc(t), rt(m), minfr(m), maxfr(m), frs(m), CV, fs );
    end
    fprintf('%d Quasistatic inverse model: spike trains generated\n', m);
end
clear prev_state t m


%% IPI signal generation out of spikes signal (for gain nonlinearity)
[~, lq_ipi] = sawtooth2ipi(spikes2sawtooth([lq_spikes(2:end, :); zeros(1,N)]));
% Identification data
lq_gain = ones(size(lq_spikes));
for i = 1:N
    lq_gain(:,i) = get_force_gain_fuglevand(lq_ipi(:,i), Tf(i));
end


%% Generate identification and validation force
lq_response = zeros(T,1); % Identification data
for i = 1:N
    for t = 1:T % - max_twitch_len - 1
        if lq_spikes(t,i)
            twitch_to_add = twitches_cell{i};
            to_take = min(length(twitch_to_add), T - t);
            lq_response(t+1:(t+to_take))  = ...
                    lq_response(t+1:(t+to_take)) + (1/fmax) * lq_gain(t,i) * twitch_to_add(1:to_take);
        end
    end
    fprintf('%d Quasistatic inverse model: twitch trains are generated\n', i);
end


%% Invert the curve, get polynomial approximate
lqi_model = polyfit(lq_response, lq_exc, 3);

figure; plot(lq_response, lq_exc); hold on; plot(lq_response, polyval(lqi_model, lq_response), 'linewidth', 2);
axis tight
ylabel('Excitation');
xlabel('Force, normalized');
title('Force response to quasistatic linearly increasing excitation and its polynomial fit');
legend('Quasistatic response', 'Polynomial fit', 'Location', 'southeast'); 