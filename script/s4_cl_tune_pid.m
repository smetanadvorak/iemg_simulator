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
sys_ident_spikes = mu_pool.mn_pool.generate_spike_train_gauss(1:sys_T, nan(mu_pool.N,1), sys_ident_excitation, fs);
sys_valid_spikes = mu_pool.mn_pool.generate_spike_train_gauss(1:sys_T, nan(mu_pool.N,1), sys_valid_excitation, fs);


%% Generate identification and validation force
sys_ident_force = mf_mdl.generate_force_offline(sys_ident_spikes);
sys_valid_force = mf_mdl.generate_force_offline(sys_valid_spikes);


%% Show identification and validation data
figure; 
subplot(2,1,1); plot(sys_timeline,sys_ident_force); hold on; 
plot(sys_timeline,sys_ident_excitation); xlabel('Time, s');
title('Identification data'); legend('Force', 'Excitation'); axis tight;
subplot(2,1,2); plot(sys_timeline,sys_valid_force); hold on; 
plot(sys_timeline,sys_valid_excitation); xlabel('Time, s');
title('Validation data'); legend('Force', 'Excitation'); axis tight;


%% Set up the system identification data
sys_ident_data = iddata(sys_ident_force, sys_ident_excitation, 1/fs); % Create an identification data object
sys_ident_data.InputName = 'Neural excitation';
sys_ident_data.OutputName = 'Force';
sys_ident_data.TimeUnit = 'Seconds';
sys_ident_data.InputUnit = 'Normalized';
sys_ident_data.OutputUnit = 'Normalized';

sys_valid_data = iddata(sys_valid_force, sys_valid_excitation, 1/fs); % Create an identification data object
sys_valid_data.InputName = 'Neural excitation';
sys_valid_data.OutputName = 'Force';
sys_valid_data.TimeUnit = 'Seconds';
sys_valid_data.InputUnit = 'Normalized';
sys_valid_data.OutputUnit = 'Normalized';

%% Downsample the iddata
sys_ident_data_r = resample(sys_ident_data, 1, round(fs/fsl));
sys_valid_data_r = resample(sys_valid_data, 1, round(fs/fsl));

%% Estimate the plant model

%e2fModel = arx(detrend(pid_ident_data),[2 1 75/1000 * fs]);
%e2fModel = oe(detrend(pid_ident_data),[1 2 75/1000 * fs]);
e2fModel = oe(sys_ident_data_r, [1 1 0]);
pole(e2fModel)
abs(pole(e2fModel))

%% Compare the estimated model's output with the estimation data
figure; 
%compare(detrend(pid_valid_data), e2fModel);
plot(sys_valid_data_r.SamplingInstants, sys_ident_data_r.OutputData); hold all;
pid_pred = predict(e2fModel, sys_ident_data_r, 0);
plot(sys_valid_data_r.SamplingInstants, pid_pred.y, '--r');
xlim([-inf, inf]); ylim([0, 1.2]); xlabel('Time, s'); ylabel('Force, normalized');
legend('Validation ouput', 'Estimated model output'); 

%% Set up and tune the PID controller
%e2fModel = upsample(tf(e2fModel),fsl);
pidc = pidtune(e2fModel, 'PI', 8);

%% Test the PID on the linear model

% Define the closed-loop transfer function for the PID and plant model 
pid_test = feedback(pidc*e2fModel, 1);

% Draw step response (maximal excitation)
figure; step(pid_test);


%% Clear and junk
%clear pid_* sys_*


