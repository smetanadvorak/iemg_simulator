%% Generate trapezoidal profile
profile_len = 25; % simulation time seconds;
silence = 1;
percent_mvc = 100;%
slope = 10; % percent per second

profile = ContractionProfile(profile_len, fs, 'trapezoidal', percent_mvc/100, slope/100, silence);
%profile = ContractionProfile(profile_len, fs, 'constant', percent_mvc/100);
%profile = ContractionProfile(profile_len, fs, 'ramp', percent_mvc/100, slope/100, silence);
%profile = ContractionProfile(profile_len, fs, 'rect', percent_mvc/100);


%% Show generated profile
profile.show();
clear profile_len percent_mvc slope


%% Electrode shifts
traj_parameter_map = (1:profile.T)'/profile.T; % Time-dependent translation of the electrode
traj_parameter_map_type = 'Time-dependent';

%traj_parameter_map = zeros(profile.T,1);
%traj_parameter_map_type = 'No motion';


