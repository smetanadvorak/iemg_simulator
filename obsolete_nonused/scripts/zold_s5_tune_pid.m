T = 32*fs;

input_levels = linspace(0, 1, 1); %
input_period = fs/2;
warning('off','Ident:dataprocess:idinput7');
pseudorandom_input = idinput(2*T*numel(input_levels), 'prbs',[0 1/input_period],[0,1]);
warning('on','Ident:dataprocess:idinput7');

pseudorandom_input = sum(reshape(pseudorandom_input, [], numel(input_levels)),2); % Generate multilevel signal
pseudorandom_input = majority_filter(pseudorandom_input, 2*round(input_period/8)+1)/numel(input_levels);

l_input = round(length(pseudorandom_input)/2);
sys_ident_excitation = pseudorandom_input(1:l_input);
sys_valid_excitation = pseudorandom_input(l_input+1:end);

sys_ident_excitation(1:input_period) = 0;
sys_valid_excitation(1:input_period) = 0;
timeline = linspace(0,numel(sys_ident_excitation)/fs, numel(sys_ident_excitation));
%clear profile_nodes profile_vals pseudorandom_input


%% Generate identification and validation input
pid_ident_spikes = zeros(T,N); % Identification data
for m = 1:N
    ident_prev_state = nan;
    for t = 1:T
        [pid_ident_spikes(t,m), ident_prev_state] = generate_spike_train( t, ident_prev_state, sys_ident_excitation(t), rt(m), minfr(m), maxfr(m), frs(m), CV, fs );
    end
    fprintf('%d Identification spike trains generated\n', m);
end

pid_valid_spikes = zeros(T,N); % Validation data
for m = 1:N
    valid_prev_state = nan;
    for t = 1:T
        [pid_valid_spikes(t,m), valid_prev_state] = generate_spike_train( t, valid_prev_state, sys_valid_excitation(t), rt(m), minfr(m), maxfr(m), frs(m), CV, fs );
    end
    fprintf('%d Validation spike trains generated\n', m);
end
clear valid_prev_state


%% IPI signal generation out of spikes signal (for gain nonlinearity)
[~, pid_ident_ipi] = sawtooth2ipi(spikes2sawtooth([pid_ident_spikes(2:end, :); zeros(1,N)]));

% Identification data
pid_ident_gain = ones(size(pid_ident_spikes));
for i = 1:N
    pid_ident_gain(:,i) = get_force_gain_fuglevand(pid_ident_ipi(:,i), Tf(i));
end

% Validation data
[~, pid_valid_ipi] = sawtooth2ipi(spikes2sawtooth([pid_valid_spikes(2:end, :); zeros(1,N)]));
pid_valid_gain = ones(size(pid_valid_spikes));
for i = 1:N
    pid_valid_gain(:,i) = get_force_gain_fuglevand(pid_valid_ipi(:,i), Tf(i));
end


%% Generate identification and validation force
pid_ident_force = zeros(T, 1); % Identification data
for i = 1:N
    for t = 1:T% - max_twitch_len - 1
        if pid_ident_spikes(t,i)
            twitch_to_add = twitches_cell{i};
            to_take = min(length(twitch_to_add), T - t);
            pid_ident_force(t:(t+to_take-1))  = ...
                    pid_ident_force(t:(t+to_take-1)) + (1/fmax) * pid_ident_gain(t,i) * twitch_to_add(1:to_take);
        end
    end
    fprintf('%d Identification twitch trains are generated\n', i);
end

pid_valid_force = zeros(T, 1); % Validation data
for i = 1:N
    for t = 1:T% - max_twitch_len - 1
        if pid_valid_spikes(t,i)
            twitch_to_add = twitches_cell{i};
            to_take = min(length(twitch_to_add), T - t);
            pid_valid_force(t:(t+to_take-1))  = ...
                    pid_valid_force(t:(t+to_take-1)) + (1/fmax) * pid_valid_gain(t,i) * twitch_to_add(1:to_take);
        end
    end
    fprintf('%d Validation twitch trains are generated\n', i);
end
clear to_take twitch_to_add i t



%% Show identification and validation data
figure; 
subplot(2,1,1); plot(timeline,pid_ident_force); hold on; 
plot(timeline,sys_ident_excitation); xlabel('Time, s');
title('Identification data'); legend('Force', 'Excitation'); axis tight;
subplot(2,1,2); plot(timeline,pid_valid_force); hold on; 
plot(timeline,sys_valid_excitation); xlabel('Time, s');
title('Validation data'); legend('Force', 'Excitation'); axis tight;



%% Set up the system identification data
pid_ident_data = iddata(pid_ident_force, sys_ident_excitation); % Create an identification data object
pid_ident_data.InputName = 'Neural excitation';
pid_ident_data.OutputName = 'Force';
pid_ident_data.TimeUnit = 'Seconds';
pid_ident_data.InputUnit = 'Normalized';
pid_ident_data.OutputUnit = 'Normalized';

pid_valid_data = iddata(pid_valid_force, sys_valid_excitation); % Create an identification data object
pid_valid_data.InputName = 'Neural excitation';
pid_valid_data.OutputName = 'Force';
pid_valid_data.TimeUnit = 'Seconds';
pid_valid_data.InputUnit = 'Normalized';
pid_valid_data.OutputUnit = 'Normalized';



%% Estimate the plant model

%e2fModel = arx(detrend(pid_ident_data),[2 1 75/1000 * fs]);
%e2fModel = oe(detrend(pid_ident_data),[1 2 75/1000 * fs]);
e2fModel = oe(pid_ident_data, [1 2 0]);

%% Compare the estimated model's output with the estimation data
figure; 
%compare(detrend(pid_valid_data), e2fModel);
plot(pid_ident_data.SamplingInstants, pid_ident_data.OutputData); hold all;
pid_pred = predict(e2fModel, pid_ident_data, 0);
plot(pid_pred, '--r');
legend('Identification ouput', 'Estimated model output');

%% Set up and tune the PID controller
pidc = pidtune(e2fModel, 'PI');



%% Test the PID on the linear model

% Define the closed-loop transfer function for the PID and plant model 
pid_test = feedback(pidc*e2fModel, 1);

% Draw step response (maximal excitation)
figure; step(pid_test);


%% Clear and junk
%clear pid_*


