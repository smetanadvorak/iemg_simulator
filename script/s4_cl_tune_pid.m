sys_T = 32*fs;

sys_pseudorandom_input = multilevel_prbs(sys_T, fs/2, 1);

sys_l_input = round(length(sys_pseudorandom_input)/2);
sys_ident_excitation = sys_pseudorandom_input(1:sys_l_input);
sys_valid_excitation = sys_pseudorandom_input(sys_l_input+1:end);

sys_ident_excitation(1:round(fs/2)) = 0;
sys_valid_excitation(1:round(fs/2)) = 0;
sys_timeline = linspace(0,numel(sys_ident_excitation)/fs, numel(sys_ident_excitation));

clear l_input 
%% Generate identification and validation input
pid_ident_spikes = mu_pool.mn_pool.generate_spike_train_gauss(1:sys_T, nan(N,1), sys_ident_excitation, fs);
pid_valid_spikes = mu_pool.mn_pool.generate_spike_train_gauss(1:sys_T, nan(N,1), sys_valid_excitation, fs);


%% Generate identification and validation force
pid_ident_force = mf_mdl.generate_force_offline(pid_ident_spikes);
pid_valid_force = mf_mdl.generate_force_offline(pid_valid_spikes);


%% Show identification and validation data
figure; 
subplot(2,1,1); plot(sys_timeline,pid_ident_force); hold on; 
plot(sys_timeline,sys_ident_excitation); xlabel('Time, s');
title('Identification data'); legend('Force', 'Excitation'); axis tight;
subplot(2,1,2); plot(sys_timeline,pid_valid_force); hold on; 
plot(sys_timeline,sys_valid_excitation); xlabel('Time, s');
title('Validation data'); legend('Force', 'Excitation'); axis tight;


%% Set up the system identification data
pid_ident_data = iddata(pid_ident_force, sys_ident_excitation, 1/fs); % Create an identification data object
pid_ident_data.InputName = 'Neural excitation';
pid_ident_data.OutputName = 'Force';
pid_ident_data.TimeUnit = 'Seconds';
pid_ident_data.InputUnit = 'Normalized';
pid_ident_data.OutputUnit = 'Normalized';

pid_valid_data = iddata(pid_valid_force, sys_valid_excitation, 1/fs); % Create an identification data object
pid_valid_data.InputName = 'Neural excitation';
pid_valid_data.OutputName = 'Force';
pid_valid_data.TimeUnit = 'Seconds';
pid_valid_data.InputUnit = 'Normalized';
pid_valid_data.OutputUnit = 'Normalized';

%% Downsample the iddata
pid_ident_data_r = resample(pid_ident_data, 1, round(fs/fsl));
pid_valid_data_r = resample(pid_valid_data, 1, round(fs/fsl));

%% Estimate the plant model

%e2fModel = arx(detrend(pid_ident_data),[2 1 75/1000 * fs]);
%e2fModel = oe(detrend(pid_ident_data),[1 2 75/1000 * fs]);
e2fModel = oe(pid_ident_data_r, [1 2 0]);
pole(e2fModel)
abs(pole(e2fModel))

%% Compare the estimated model's output with the estimation data
figure; 
%compare(detrend(pid_valid_data), e2fModel);
plot(pid_ident_data_r.SamplingInstants, pid_ident_data_r.OutputData); hold all;
pid_pred = predict(e2fModel, pid_ident_data_r, 0);
plot(pid_pred, '--r');
legend('Identification ouput', 'Estimated model output');

%% Set up and tune the PID controller
%e2fModel = upsample(tf(e2fModel),fsl);
pidc = pidtune(e2fModel, 'PI', 8);

%% Test the PID on the linear model

% Define the closed-loop transfer function for the PID and plant model 
pid_test = feedback(pidc*e2fModel, 1);

% Draw step response (maximal excitation)
figure; step(pid_test);


%% Clear and junk
clear pid_* sys_*


