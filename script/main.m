
run init    % If cd(fileparts(mfilename('fullpath'))) fail to execute, choose the
            % script/ folder as current matlab folder and run the rest of init.m

fs = 10000; % Time sampling frequency for EMG, [Hz]
dt = 1/fs; % Sampling period of EMG [s]
fsl = 50; % Time sampling frequency for force, [Hz]
dz = 0.5; % Spatial sampling frequency for muscle fibers action potentials [mm]   

%% MU pool model
run s01_cl_init_mn_pool  % Generate the motor neuron pool (defines the muscle geometry,
                        % as well as sizes and innervation centers of the motor neurons)
run s02_cl_init_mu_pool  % Generate the muscle fibers coordinates, assign the muscle fibers to
                        % the MNs to generate motor units. Assign diameters of fibers. 

%% Muscle model
run s03_cl_init_force_model  % Generate twitch waveforms and non-linearity coefficient (follows Fuglevand 1993)
run s04_cl_tune_pid          % Identifies the excitation-force model. Establishes a PID controller that outputs the
                            % excitation needed to follow the force profile.

%% Simulation of MUAPs (may be specifically time-consuming)
run s07_cl_init_electrode    % Defines the electrode's geometry (form and number of channels) and position in the muscle
run s08_cl_init_muaps        % Defines the terminal arborization geometry and pre-calculates the MUAPs
run s09_cl_generate_mvc_emg % Generates the EMG at maximum voluntary contraction, for the noise level reference
run s17_cl_save_model       % Saves the model to a .m file


%% Simulation of contraction: force and EMG generation
run s05_cl_init_profile      % Define the contraction profile (trapezoidal, constant, etc.)
run s06_cl_generate_force    % Generate the spike trains according to force profile
run s10_cl_generate_emg     % Generate the EMG out of the spike trains.


%% Annotation for decomposition
run s11_cl_get_detectable_mus   % Decides which MUAPs are going to the dictionary.
run s12_cl_generate_annotation  % Generates the full annotation (spike trains of all MN)
run s13_cl_generate_dictionary  % Generates the MUAPs dictionary
run s14_cl_generate_annotation_for_decomp  % Generate the decomposition annotation (spike trains of MNs in dictionary)
run s15_cl_reconstruct_signal   % Reconstructs EMG from decomposition annotation and dictionary, to check if everything is fine.

%% Save the generated model
run s18_cl_save_output          % Saves the output signals and annotations to the disk
run s16_cl_make_a_report        % Generates a .txt file with model description


