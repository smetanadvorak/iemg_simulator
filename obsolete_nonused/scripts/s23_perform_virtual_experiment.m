spikes = zeros(T,N);
decomposition = zeros(T,numel(detectable_ind));
force_estimate = zeros(T,1);
excitation = zeros(T,1);
estimator.reset();

err = zeros(T,1);
int_err = 0;
der_err = 0;

Kp = pidc.Kp;
Ki = pidc.Ki * pidc.Ts;
Kd = pidc.Kd;

prev_neuron_state = nan(N,1);
for t = 1:T-1
    % Error calculation for PID controller's input
    err(t) = virtual_profile(t) - force_estimate(t); % Proportional error
    if virtual_profile(t) > 0, int_err = int_err + err(t); else int_err = 0; end % Integrated error (supressed when no goal force)
    if t > 1, der_err = err(t) - err(t-1); end % Local derivative of the error 
    
    % PID controller's output is net excitation to the motor neuron pool
	excitation(t) = Kp * err(t) + Ki * int_err + Kd * der_err;
    excitation(t) = max(excitation(t), 0);
    excitation(t) = min(excitation(t), 1);
    
    % Spike train generation
    for m = 1:N
        [spikes(t,m), prev_neuron_state(m), ipi] = generate_spike_train(t, prev_neuron_state(m), excitation(t), rt(m), minfr(m), maxfr(m), frs(m), CV, fs);
    end
    
    % Decomposition generation
    decomposition(t, :) = spikes(t, detectable_ind);
    
    % Force estimate calculation
    estimator.update(decomposition(t,:));
    force_estimate(t+1) = estimator.predict();
    
    if ~mod(t,fs)
        fprintf('%d seconds generated\n', floor(t/fs));
    end
end
clear t m twitch_to_add to_take gain int_err der_err


%% Plot goal force + resulting force, excitation and error
figure; hold all; 
plot(profile_timeline, virtual_profile, 'k--', 'linewidth', 2); 
plot(profile_timeline, force_estimate,'g'); 
plot(profile_timeline, err, 'r');
ylabel('Force, normalized');
%yyaxis right;
%plot(profile_timeline, excitation,'b'); 
%ylabel('Excitation');
legend('Goal force', 'Force estimate', 'Error', 'Resulting excitation');


%% Plot resulting spikes
firings = spikes2firings(spikes);
figure;
for m = 1:length(firings)
    indices = firings{m}/fs;
    separator = repmat(m, length(indices), 1);
    if ismember(m, detectable_ind)
        plot(indices, separator, 'go', 'markersize',3); hold on;
    else
        plot(indices, separator, 'ko', 'markersize',3); hold on;
    end
end

ylabel('Motor neuron index');
yyaxis right;
plot(profile_timeline, excitation, 'linewidth', 2);
ylabel('Excitation, normalized units');
xlabel('Time, sec');
title('Motor neurons'' firining times and excitation vs time');

