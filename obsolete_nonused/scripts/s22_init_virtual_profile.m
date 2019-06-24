%% Generate trapezoidal profile
T = 10 * fs; % simulation time in samples;
profile_timeline = (1:T)/fs;

percent_mec = 25; %

%% Trapezoidal profile
profile_nodes = [0, 5, 40, 60, 95, 100]/100 * T/fs; % trapezoidal profile
profile_vals = [0, 0, 1 + eps, 1 + eps, 0, 0]; % trapezoidal profile

profile_vals = profile_vals * 1 * MEC * percent_mec/100; %+eps to trigger the last motor neuron
virtual_profile = interp1(profile_nodes, profile_vals, profile_timeline);
profile_type = ['trapezoidal_' num2str(percent_mec) '%MVC'];

%% Constant profile
% virtual_profile = (percent_mec/100) * MEC *ones(size(profile_timeline));
% profile_type = ['constant_', num2str(percent_mvc), '%MVC'];

%% Ramp profile
% profile_nodes = [0, 10, 90, 100]/100 * T/fs; % trapezoidal profile
% profile_vals = [0, 0, 1 + eps, 1 + eps]; % trapezoidal profile
% 
% profile_vals = profile_vals * MEC * percent_mvc/100; %+eps to trigger the last motor neuron
% virtual_profile = interp1(profile_nodes, profile_vals, profile_timeline);
% profile_type = ['ramp_' num2str(percent_mvc) '%MVC'];


%% Electrode shifts
% traj_parameter_map = (1:T)'/T; % Time-dependent translation of the electrode
% traj_parameter_map_type = 'Time-dependent';

traj_parameter_map = zeros(T,1);
traj_parameter_map_type = 'No motion';


clear profile_vals profile_nodes
