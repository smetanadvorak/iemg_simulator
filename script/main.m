cd(fileparts(mfilename('fullpath')));
addpath(pwd);
work_folder = '../';
cd(work_folder);
addpath(genpath('./lib'));
addpath(genpath('./dependancies'));

fs = 20000; % Time sampling frequency for EMG, [Hz]
fsl = 50; % Time sampling frequency for force
dt = 1/fs; % [s]
dz = 0.5; % [mm]   

%% MU pool
run s1_cl_init_mn_pool % Attention, deluca model implemented. Check outputs next time you simulate
run s2_cl_init_mu_pool 

%% Muscle
run s3_cl_init_force_model
run s4_cl_tune_pid

%% MUAPs and EMG (long)
run s8_cl_init_electrode
run s9_cl_init_muaps
run s10_cl_generate_mvc_emg
run s17_cl_save_model

%% Contraction
% run s5_cl_init_profile
% run s6_cl_generate_force
% run s11_cl_generate_emg

%% Annotation
% run s12_cl_get_detectable_mus
% run s12_cl_generate_annotation
% run s13_cl_generate_dictionary
% run s14_cl_generate_annotation_for_decomp
% run s14_cl_reconstruct_signal

%% Save the generated model
% run s18_cl_save_output
% run s15_cl_make_a_report


