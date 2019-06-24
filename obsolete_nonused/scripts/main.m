%% Integration parameters
% Attentions:
% 1) rad_hood is switched off in assign_mf_gauss_size.m
% 2) distances from observation points to fibers are lower bounded by mean
% diameter (kind of has sense).
% 3) 

% ToDos:
% 1) Dictionaries output (at timestamps?).
% 2) 
addpath(pwd);
work_folder = '../';
cd(work_folder);
addpath(genpath('./lib'));
addpath(genpath('./dependancies'));

fs = 25000; % EMG sampling frequency, [Hz]
fsl = 25;   % Mechanics sampling frequency [Hz]
dt = 1/fs; % [s]
dz = 0.5; % [mm] 

CV = 1/6;
    
%% Init necessary
run s1_init_pool
run s2_init_muscle

%% Excitation-force controller
run s3_init_force_model
run s5_tune_pid

%% Add EMG simulation
run s8_init_electrode
run s9_init_muaps
run s10_generate_mvc_emg

%%
run s6_init_profile
run s7_generate_force

%%
run s11_generate_emg
run s12_get_detectable_mus
run s12_generate_annotation
run s13_generate_dictionary
run s14_reconstruct_signal
run s16_save_output
run s15_make_a_report



%%
run s17_save_model
