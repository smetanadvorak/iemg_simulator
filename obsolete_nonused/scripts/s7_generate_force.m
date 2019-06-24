spikes = zeros(T,N);
ipi = zeros(T,N);
force = zeros(T,1);
excitation = zeros(T,1);

error = zeros(T,1);
int_err = 0;
der_err = 0;

Kp = pidc.Kp;
Ki = pidc.Ki * pidc.Ts;
Kd = pidc.Kd;

prev_neuron_state = nan(N,1);

sub_counter = 0;

sub_profile = force_profile((1 : round(T/fs*fsl)) * round(fs/fsl));
sub_error = zeros(round(T/fs*fsl), 1);
sub_force = zeros(round(T/fs*fsl), 1);
sub_exc = zeros(round(T/fs*fsl), 1);
sub_t = 1;

mf_mdl.init_online_buffer();
for t = 1:T-1
    
    % Error accumulation for subsampling
    error(t) = force_profile(t) - force(t);
    
    sub_counter = sub_counter + 1;
    
    if ~mod(sub_counter, round(fs/fsl))     
        % Error calculation for PID controller's input
        sub_force(sub_t) = mean(force(t-fsl+1:t));
        sub_error(sub_t) = sub_profile(sub_t) - sub_force(sub_t); % Proportional error
        if sub_profile(sub_t) > 0, int_err = int_err + sub_error(sub_t); else int_err = 0; end % Integrated error (supressed when no goal force)
        if sub_t > 1, der_err = sub_error(sub_t) - sub_error(sub_t-1); end % Local derivative of the error
        
        % PID controller's output is net excitation to the motor neuron pool
        sub_exc(sub_t) = Kp * sub_error(sub_t) + Ki * int_err + Kd * der_err;
        sub_exc(sub_t) = max(sub_exc(sub_t), 0);
        sub_exc(sub_t) = min(sub_exc(sub_t), 1);
        
        % The inverse model
        if sub_profile(sub_t) <= 0
            sub_exc(sub_t) = 0;
        elseif sub_profile(sub_t) >= 1
            sub_exc(sub_t) = 1;
        else
            sub_exc(sub_t) = sub_exc(sub_t) + polyval(qsi_model, sub_profile(sub_t));
        end
        
        sub_t = sub_t + 1;
    end
    
    % Upsampling back...
    if sub_t == 1
        excitation(t) = 0;
    else
        excitation(t) = sub_exc(sub_t-1);
    end
    
    % Force generation
    [spikes(t,:), prev_neuron_state, ipi(t,:)] = generate_spike_train_vect(t, prev_neuron_state, excitation(t), rt, minfr, maxfr, frs, CV, fs);
    force(t+1) = mf_mdl.generate_force_online(spikes(t,:), ipi(t,:));
    
    if ~mod(t,fs)
        fprintf('%d seconds generated\n', floor(t/fs));
    end
end

%% Plot goal force + resulting force, excitation and error
figure; hold all; 
plot(profile_timeline, force_profile * 100, 'k'); 
plot(profile_timeline, force * 100,'g'); 
plot(profile_timeline, error * 100, 'r');

sub_timeline = (1 : round(T/fs*fsl)) * round(fs/fsl) / fs;
plot(sub_timeline, sub_error * 100, 'r.', 'linestyle', 'none');
plot(sub_timeline, sub_force * 100, 'g.', 'linestyle', 'none');

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

clear t m twitch_to_add to_take gain int_err der_err sub_* Kd Kp Ki int_err der_err prev_neuron_state