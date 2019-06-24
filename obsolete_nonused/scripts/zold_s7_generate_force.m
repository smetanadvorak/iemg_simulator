spikes = zeros(T,N);
ipi = zeros(T,N);
force = zeros(T,1);
excitation = zeros(T,1);

err = zeros(T,1);
int_err = 0;
der_err = 0;

Kp = pidc.Kp;
Ki = pidc.Ki * pidc.Ts;
Kd = pidc.Kd;

prev_neuron_state = nan(N,1);
for t = 1:T
    % Error calculation for PID controller's input
    err(t) = force_profile(t) - force(t); % Proportional error
    if force_profile(t) > 0, int_err = int_err + err(t); else int_err = 0; end; % Integrated error (supressed when no goal force)
    if t > 1, der_err = err(t) - err(t-1); end; % Local derivative of the error 
    
    % PID controller's output is net excitation to the motor neuron pool
	excitation(t) = Kp * err(t) + Ki * int_err + Kd * der_err;
    excitation(t) = max(excitation(t), 0);
    excitation(t) = min(excitation(t), 1);
    
    excitation(t) = excitation(t) + polyval(qsi_model, force_profile(t));
    
    % Force generation
    for m = 1:N
        [spikes(t,m), prev_neuron_state(m), ipi(t,m)] = generate_spike_train(t, prev_neuron_state(m), excitation(t), rt(m), minfr(m), maxfr(m), frs(m), CV, fs);
        if spikes(t,m)
            twitch_to_add = twitches_cell{m};
            to_take = min(length(twitch_to_add), T - t);
            gain = get_force_gain_fuglevand(ipi(t,m), Tf(m));
            force(t:(t+to_take-1)) = force(t:(t+to_take-1)) + (1/fmax) * gain * twitch_to_add(1:to_take);
        end
    end
    if ~mod(t,fs)
        fprintf('%d seconds generated\n', floor(t/fs));
    end
end
clear t m twitch_to_add to_take gain int_err der_err


%% Plot goal force + resulting force, excitation and error
figure; hold all; 
plot(profile_timeline, force_profile * 100, 'k'); 
plot(profile_timeline, force * 100,'g'); 
plot(profile_timeline, err * 100, 'r');
ylabel('Force, \%MVC');
xlabel('Time, s');
%yyaxis right;
%plot(profile_timeline, excitation,'b'); 
%ylabel('Excitation');
legend('Goal force', 'Resulting force', 'Error');%, 'Resulting excitation');


%% Plot resulting spikes
firings = spikes2firings(spikes);
figure;
for m = 1:length(firings)
    indices = firings{m}/fs;
    separator = repmat(m, length(indices), 1);
    plot(indices, separator, 'ko'); hold on;
end
ylabel('Motor neuron index');
yyaxis right;
plot(profile_timeline, excitation, 'linewidth', 2);
ylabel('Excitation, normalized units');
xlabel('Time, sec');
title('Motor neurons'' firining times and excitation vs time');

