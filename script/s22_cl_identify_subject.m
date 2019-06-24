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
ident_input = pseudorandom_input(1:l_input);
valid_input = pseudorandom_input(l_input+1:end);

ident_input(1:input_period) = 0;
valid_input(1:input_period) = 0;
timeline = linspace(0,numel(ident_input)/fs, numel(ident_input));

% Scale to the maximal excitation
ident_input = ident_input * max_excitation;
valid_input = valid_input * max_excitation;

% Calculate spike trains
ident_spikes = mu_pool.mn_pool.generate_spike_train_gauss(1:T_ident, nan(mu_pool.N,1), ident_input, fs);
valid_spikes = mu_pool.mn_pool.generate_spike_train_gauss(1:T_ident, nan(mu_pool.N,1), valid_input, fs);

ident_spikes = ident_spikes(:,detectable_ind);
valid_spikes = valid_spikes(:,detectable_ind);

%% Downsample the ident data
T_ident = T_ident/fs*fsl;

downsampler_aux.reset();
downsampler_spk.reset();
downsampler_act.reset();

ident_input_dwns = zeros(size(ident_input));
ident_act_dwns = zeros(size(ident_spikes));
ident_spk_dwns  = zeros(size(ident_spikes));
ident_act = spikes2activity(ident_spikes, length(fe_window)/fsl*fs);

for t = 1:size(ident_spikes,1)
    ident_input_dwns(t) = downsampler_aux.update(ident_input(t));
    ident_act_dwns(t,:) = downsampler_act.update(ident_act(t,:));
    ident_spk_dwns(t,:) = downsampler_spk.update(ident_spikes(t,:));
end

ident_input_dwns(isnan(ident_input_dwns(:,1)),:) = [];
ident_act_dwns(isnan(ident_act_dwns(:,1)),:) = [];
ident_spk_dwns(isnan(ident_spk_dwns(:,1)),:) = [];

ident_act_dwns = double(ident_act_dwns >= 0.5);
ident_spk_dwns = double(ident_spk_dwns > 0);

% Calculate firing rates
ident_fr_dwns = zeros(size(ident_spk_dwns));
for m = 1:size(ident_spk_dwns,2)
    fr_estimator(m).reset();
    for t = 1:size(ident_spk_dwns,1)
        ident_fr_dwns(t,m) = fr_estimator(m).update(ident_spk_dwns(t,m)) * fsl;
    end
    fr_estimator(m).reset();
end

T_ident = size(ident_input_dwns,1);

%% Calculate the output for the ident data
fprintf('Calculating the estimates for identification and validation data');
sys_ident_output = zeros(T_ident, 1); % Identification data
estimator.reset();

switch estimator.type
    case 'RS'
        for i = 1:numel(sys_ident_output)
            estimator.update(ident_act_dwns(i,:));
            sys_ident_output(i) = estimator.predict();
        end
        
    case 'RC'
        for i = 1:numel(sys_ident_output)
            estimator.update(ident_fr_dwns(i,:));
            sys_ident_output(i) = estimator.predict();
        end
        
    case 'RSRC'
        for i = 1:numel(sys_ident_output)
            estimator.update(ident_fr_dwns, ident_act_dwns(i,:));
            sys_ident_output(i) = estimator.predict();
        end
        
    case 'CDR'
        for i = 1:numel(sys_ident_output)
            estimator.update(ident_spk_dwns(i,:));
            sys_ident_output(i) = estimator.predict();
        end
end

%% Downsample the validation data
downsampler_aux.reset();
downsampler_spk.reset();
downsampler_act.reset();

valid_input_dwns = zeros(size(valid_input));
valid_act_dwns = zeros(size(valid_spikes));
valid_spk_dwns  = zeros(size(valid_spikes));
valid_act = spikes2activity(valid_spikes, length(fe_window)/fsl*fs);

for t = 1:size(valid_spikes,1)
    valid_input_dwns(t) = downsampler_aux.update(valid_input(t));
    valid_act_dwns(t,:) = downsampler_act.update(valid_act(t,:));
    valid_spk_dwns(t,:) = downsampler_spk.update(valid_spikes(t,:));
end

valid_input_dwns(isnan(valid_input_dwns(:,1)),:) = [];
valid_act_dwns(isnan(valid_act_dwns(:,1)),:) = [];
valid_spk_dwns(isnan(valid_spk_dwns(:,1)),:) = [];

valid_act_dwns = double(valid_act_dwns >= 0.5);
valid_spk_dwns = double(valid_spk_dwns > 0);

% Calculate firing rates
valid_fr_dwns = zeros(size(valid_spk_dwns));
for m = 1:size(valid_spk_dwns,2)
    fr_estimator(m).reset();
    for t = 1:size(valid_spk_dwns,1)
        valid_fr_dwns(t,m) = fr_estimator(m).update(valid_spk_dwns(t,m)) * fsl;
    end
    fr_estimator(m).reset();
end

T_ident = size(valid_input_dwns,1);

%% Calculate the output of the simulator + estimator
fprintf('Calculating the estimates for identification and validation data');
sys_valid_output = zeros(T_ident, 1); % Identification data
estimator.reset();

switch estimator.type
    case 'RS'
        for i = 1:numel(sys_valid_output)
            estimator.update(valid_act_dwns(i,:));
            sys_valid_output(i) = estimator.predict();
        end
        
    case 'RC'
        for i = 1:numel(sys_valid_output)
            estimator.update(valid_fr_dwns(i,:));
            sys_valid_output(i) = estimator.predict();
        end
        
    case 'RSRC'
        for i = 1:numel(sys_valid_output)
            estimator.update(valid_fr_dwns, valid_act_dwns(i,:));
            sys_valid_output(i) = estimator.predict();
        end
        
    case 'CDR'
        for i = 1:numel(sys_valid_output)
            estimator.update(valid_spk_dwns(i,:));
            sys_valid_output(i) = estimator.predict();
        end
end

%% Show identification and validation data
figure; 
subplot(2,1,1); plot(sys_ident_output); hold on; plot(ident_input_dwns); title('Identification data'); legend('Force', 'Excitation');
subplot(2,1,2); plot(sys_valid_output); hold on; plot(valid_input_dwns); title('Validation data'); legend('Force', 'Excitation');



%% Estimate the system
sys_ident_data = iddata(sys_ident_output, ident_input_dwns, 1/fs); % Create an identification data object
sys_ident_data.InputName = 'Neural excitation';
sys_ident_data.OutputName = 'Force';
sys_ident_data.TimeUnit = 'Seconds';
sys_ident_data.InputUnit = 'Normalized';
sys_ident_data.OutputUnit = 'Normalized';

sys_valid_data = iddata(sys_valid_output, valid_input_dwns, 1/fs); % Create an identification data object
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



