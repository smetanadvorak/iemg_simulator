%% Set the virtual user
% Estimate the linear model of the controlled system

%% Generate the identification and validation inputs
T_ident = 32*fs;

input_levels = linspace(0, 1, 4); %
input_period = fs/2;
warning('off','Ident:dataprocess:idinput7');
pseudorandom_input = idinput(2*T_ident*numel(input_levels), 'prbs',[0 1/input_period],[0,1]);
warning('on','Ident:dataprocess:idinput7');

pseudorandom_input = sum(reshape(pseudorandom_input, [], numel(input_levels)),2); % Generate multilevel signal
pseudorandom_input = majority_filter(pseudorandom_input, 2*round(input_period/8)+1)/numel(input_levels);

l_input = round(length(pseudorandom_input)/2);
sys_ident_input = pseudorandom_input(1:l_input);
sys_valid_input = pseudorandom_input(l_input+1:end);

sys_ident_input(1:input_period) = 0;
sys_valid_input(1:input_period) = 0;
timeline = linspace(0,numel(sys_ident_input)/fs, numel(sys_ident_input));

% Scale to the maximal excitation
sys_ident_input = sys_ident_input * max_excitation;
sys_valid_input = sys_valid_input * max_excitation;
    
%sys_ident_input = generate_identification_input( 2*fs, fs, (0.2:0.2:1)*max_excitation, 0.1);
%sys_valid_input = generate_identification_input( 2*fs, fs, (0.2:0.2:1)*max_excitation, 0.1);
%T_ident = numel(sys_ident_input);

%% Calculate spike trains
sys_ident_spikes = zeros(T_ident,numel(detectable_ind)); % Identification data
for m = active_ind(:)'
    ident_prev_state = nan;
    for t = 1:T_ident
        [sys_ident_spikes(t,m), ident_prev_state] = generate_spike_train( t, ident_prev_state, sys_ident_input(t), rt(m), minfr(m), maxfr(m), frs(m), CV, fs );
    end
    fprintf('%d Identification spike trains generated\n', m);
end
%pid_ident_spikes = pid_ident_spikes(:, detectable_ind);

sys_valid_spikes = zeros(T_ident,numel(detectable_ind)); % Validation data
for m = active_ind(:)'
    valid_prev_state = nan;
    for t = 1:T_ident
        [sys_valid_spikes(t,m), valid_prev_state] = generate_spike_train( t, valid_prev_state, sys_valid_input(t), rt(m), minfr(m), maxfr(m), frs(m), CV, fs );
    end
    fprintf('%d Validation spike trains generated\n', m);
end
%pid_valid_spikes = pid_valid_spikes(:, detectable_ind);

clear valid_prev_state ident_prev_state


%% Calculate the output of the simulator + estimator
fprintf('Calculating the estimates for identification and validation data');
sys_ident_output = zeros(T_ident, 1); % Identification data
estimator.reset();
for i = 1:numel(sys_ident_output)
    estimator.update(sys_ident_spikes(i,:));
    sys_ident_output(i) = estimator.predict();
end

sys_valid_output = zeros(T_ident, 1); % Validation data
estimator.reset();
for i = 1:numel(sys_valid_output)
    estimator.update(sys_valid_spikes(i,:));
    sys_valid_output(i) = estimator.predict();
end

clear to_take twitch_to_add i t


%% Show identification and validation data
figure; 
subplot(2,1,1); plot(sys_ident_output); hold on; plot(sys_ident_input); title('Identification data'); legend('Force', 'Excitation');
subplot(2,1,2); plot(sys_valid_output); hold on; plot(sys_valid_input); title('Validation data'); legend('Force', 'Excitation');



%% Estimate the system
sys_ident_data = iddata(sys_ident_output, sys_ident_input, 1/fs); % Create an identification data object
sys_ident_data.InputName = 'Neural excitation';
sys_ident_data.OutputName = 'Force';
sys_ident_data.TimeUnit = 'Seconds';
sys_ident_data.InputUnit = 'Normalized';
sys_ident_data.OutputUnit = 'Normalized';

sys_valid_data = iddata(sys_valid_output, sys_valid_input, 1/fs); % Create an identification data object
sys_valid_data.InputName = 'Neural excitation';
sys_valid_data.OutputName = 'Force';
sys_valid_data.TimeUnit = 'Seconds';
sys_valid_data.InputUnit = 'Normalized';
sys_valid_data.OutputUnit = 'Normalized';



%% Estimate the plant model
sprintf('Identifying the model for prosthetic control ...\n');
e2fModel_virtual = oe(detrend(sys_ident_data),[1 2 0]);
figure; compare(detrend(sys_valid_data), e2fModel_virtual);
e2fModel_virtual;


%% Set up and tune the PID controller
fprintf('Tuning the PID for prosthetic control...\n');
pidc = pidtune(e2fModel_virtual, 'PI');


%% Test the PID on the linear model
fprintf('Testing the PID for prosthetic control...\n');
pid_test = feedback(pidc*e2fModel_virtual, 1);

% Draw step response (maximal excitation)
figure; step(pid_test);



% step(pid_test);
% % Create trapezoidal profile for testing
% T = 10 * fs;
% profile_timeline = (1:T)/fs;
% percent_mvc = 10; %
% profile_nodes = [0, 5, 40, 60, 95, 100]/100 * T/fs; % trapezoidal profile
% profile_vals = [0, 0, 1 + eps, 1 + eps, 0, 0]; % trapezoidal profile
% 
% profile_vals = profile_vals * percent_mvc/100; %+eps to trigger the last motor neuron
% force_profile = interp1(profile_nodes, profile_vals, profile_timeline)';
% 
% pid_valid_data = iddata(force_profile, force_profile, 1/fs);
% 
% pid_test_output = lsim(pid_test, force_profile, profile_timeline); %compare(detrend(pid_valid_data), pid_test);
% figure; plot(force_profile); hold on; plot(pid_test_output); legend('Force profile for PID test', 'Closed-loop output');

%% Clear and junk
%clear sys_*




