%% Generate trapezoidal profile
T = 10 * fs; % simulation time in samples;
profile_timeline = (1:T)/fs;

percent_mvc = 25; %

% Trapezoidal profile
profile_nodes = [0, 5, 40, 60, 95, 100]/100 * T/fs; % trapezoidal profile
profile_vals = [0, 0, 1 + eps, 1 + eps, 0, 0]; % trapezoidal profile

profile_vals = profile_vals * 1 * percent_mvc/100; %+eps to trigger the last motor neuron
force_profile = interp1(profile_nodes, profile_vals, profile_timeline);
profile_type = ['trapezoidal_' num2str(percent_mvc) '%MVC'];

%% Constant profile
% force_profile = (percent_mvc/100) * ones(size(profile_timeline));
% profile_type = ['constant_', num2str(percent_mvc), '%MVC'];

%% Ramp profile
% profile_nodes = [0, 10, 90, 100]/100 * T/fs; % trapezoidal profile
% profile_vals = [0, 0, 1 + eps, 1 + eps]; % trapezoidal profile
% 
% profile_vals = profile_vals * 1 * percent_mvc/100; %+eps to trigger the last motor neuron
% force_profile = interp1(profile_nodes, profile_vals, profile_timeline);
% profile_type = ['ramp_' num2str(percent_mvc) '%MVC'];


%% Electrode shifts
% traj_parameter_map = (1:T)'/T; % Time-dependent translation of the electrode
% traj_parameter_map_type = 'Time-dependent';

force_profile = force_profile(:);

traj_parameter_map = zeros(T,1);
traj_parameter_map_type = 'No motion';


clear profile_vals profile_nodes
