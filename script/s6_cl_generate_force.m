% Generate force according to the specified profile. Simulates the
% contraction and outputs the spike trains for this contraction.

spikes = zeros(profile.T, mu_pool.N);
ipi = zeros(profile.T,mu_pool.N);
force = zeros(profile.T,1);
excitation = zeros(profile.T,1);

error = zeros(profile.T,1);
int_err = 0;
der_err = 0;

Kp = pidc.Kp;
Ki = pidc.Ki * pidc.Ts;
Kd = pidc.Kd;

prev_neuron_state = nan(mu_pool.N,1);

sub_counter = 0;

profile.subsample(fsl);
sub_error = zeros(size(profile.sub_profile, 1), 1);
sub_force = zeros(size(profile.sub_profile,1), 1);
sub_exc = zeros(size(profile.sub_profile,1), 1);
sub_t = 1;

mf_mdl.init_online_buffer();
for t = 1:profile.T-1
    
    % Error accumulation for subsampling
    error(t) = profile.profile(t) - force(t);
    
    sub_counter = sub_counter + 1;
    
    if ~mod(sub_counter, round(fs/fsl))     
        % Error calculation for PID controller's input
        sub_force(sub_t) = mean(force(t-round(fs/fsl)+1:t));
        sub_error(sub_t) = profile.sub_profile(sub_t) - sub_force(sub_t); % Proportional error
        if profile.sub_profile(sub_t) > 0, int_err = int_err + sub_error(sub_t); else int_err = 0; end % Integrated error (supressed when no goal force)
        if sub_t > 1, der_err = sub_error(sub_t) - sub_error(sub_t-1); end % Local derivative of the error
        
        % PID controller's output is net excitation to the motor neuron pool
        sub_exc(sub_t) = Kp * sub_error(sub_t) + Ki * int_err + Kd * der_err;
        sub_exc(sub_t) = max(sub_exc(sub_t), 0);
        sub_exc(sub_t) = min(sub_exc(sub_t), 1);
        
        % The inverse model
        if profile.sub_profile(sub_t) <= 0
            sub_exc(sub_t) = 0;
        elseif profile.sub_profile(sub_t) >= 1
            sub_exc(sub_t) = 1;
        else
            sub_exc(sub_t) = sub_exc(sub_t) + mf_mdl.f2e(profile.sub_profile(sub_t));
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
    [spikes(t,:), prev_neuron_state, ipi(t,:)] = mu_pool.mn_pool.generate_spike_train_gauss(t, prev_neuron_state, excitation(t), fs);
    force(t+1) = mf_mdl.generate_force_online(spikes(t,:), ipi(t,:));
    
    if ~mod(t,fs)
        fprintf('%d seconds generated\n', floor(t/fs));
    end
end

%% Plot goal force + resulting force, excitation and error
% figure; hold all; 
% plot(profile.timeline, force * 100,'g', 'linewidth', 1.5); 
% plot(profile.timeline, error * 100, 'r', 'linewidth', 1);
% plot(profile.timeline, profile.profile * 100, 'b--', 'linewidth', 1); 
% 
% ylabel('Force, \%MVC');
% xlabel('Time, s');
% legend('Resulting force', 'Error', 'Goal force');%, 'Resulting excitation');


%% Plot resulting spikes
% firings = spikes2firings(spikes);
% figure;
% for m = 1:length(firings)
%     indices = firings{m}/fs;
%     separator = repmat(m, length(indices), 1);
%     plot(indices, separator, 'ko'); hold on;
% end
% ylabel('Motor neuron index');
% yyaxis right;
% plot(profile.timeline, excitation, 'linewidth', 2);
% ylabel('Excitation, normalized units');
% xlabel('Time, sec');
% title('Motor neurons'' firining times and excitation vs time');

%clear t m twitch_to_add to_take gain int_err der_err sub_* Kd Kp Ki int_err der_err prev_neuron_state