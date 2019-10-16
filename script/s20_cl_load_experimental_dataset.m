%% To Do:
% Make a function out of it!
% To solve: access to number of motor neurons
% Add: profile loading if generate a new dataset

%%
% This script loads the dataset for intent estimator training and
% evaluation. 

% 1) Set the dataset_path, it should contain two folders: train_data and test_data

% 2) Test/train data folders each contain a number of folders with
% contraction trials;

% 3) Each contraction trial should start and finish  with at least one 
% second of no activity (rest), in order to ensure that a trial doesn't
% affect the next trial in any way. 

% 4) Each trial folder contains at least these files:
% {anyname}.mat with EMGSIGNAL and AUXSIGNAL structures
% {anyname}.ann with the decomposition

dataset_path = '/Users/akmbpro/Nextcloud/EMG Data/ls2n/Proj_FDI_angle/FDI_Konstantin_constant_dynamic_adapted_for_estimator';
training_path = 'training_data';
testing_path = 'testing_data';

signals_file = 'simulation_signals.mat';
annotation_file = 'sim_annotation_centered_full.ann';
current_path = pwd;

try
cd(dataset_path);

%% Load training data
% Make list of trials in this folder
cd(training_path);
dir_list = dir;
dir_list = dir_list([dir_list.isdir] == 1 & ~startsWith({dir_list.name}, '.'));
dir_list = {dir_list.name};

train_emg = [];
train_aux = [];
train_spk = [];
train_max_aux = [];
fe_min_pause = 1; %In seconds;

fprintf("\nLoading training data ...");
for d = 1:numel(dir_list)
    cd(dir_list{d});
    
    % Load the signal;
    signals_file = getfield(dir('*.mat'), 'name');
    load(signals_file);
    fs = EMGSIGNAL.rate;
    fe_train_emg = EMGSIGNAL.data;
    fe_train_aux = AUXSIGNAL.data;
    fe_train_aux = filtfilt(1/round(fs/10)*ones(round(fs/10),1), 1, fe_train_aux);
    
    % Load the full annotation to extract maximal excitation (only for
    % simulated data !)
    annotation_file = getfield(dir('*.ann'), 'name');
    fe_firings = ann2firings(annotation_file, fs);
    %fe_firings{numel(fe_firings)} = [];
    fe_spikes = firings2spikes(fe_firings, length(fe_train_emg));
        
    % Append to dataset
    train_emg = [train_emg; zeros(fe_min_pause * fs, size(fe_train_emg,2)); fe_train_emg];
    train_aux = [train_aux; zeros(fe_min_pause * fs, size(fe_train_aux,2)); fe_train_aux];
    
    train_spk = [train_spk, zeros(size(train_spk,1), size(fe_spikes,2) - size(train_spk,2))];
    train_spk = [train_spk; zeros(fe_min_pause * fs, size(train_spk,2))];
    train_spk = [train_spk; fe_spikes, zeros(size(fe_spikes,1), size(train_spk,2) - size(fe_spikes,2))];
    cd ..;
end

max_active_mu = find(any(fe_spikes), 1, 'last');

MEC = max(train_aux); % Maximal estimable contraction;

cd ..

%% Load the testing data
cd(testing_path);
dir_list = dir;
dir_list = dir_list([dir_list.isdir] == 1 & ~startsWith({dir_list.name}, '.'));
dir_list = {dir_list.name};

test_emg = [];
test_aux = [];
test_spk = [];
test_max_aux = [];
fe_min_pause = 1; %In seconds;

fprintf("\nLoading testing data ...");
for d = 1:numel(dir_list)
    cd(dir_list{d});
    
    % Load the signal;
    signals_file = getfield(dir('*.mat'), 'name');
    load(signals_file);
    fs = EMGSIGNAL.rate;
    fe_test_emg = EMGSIGNAL.data;
    fe_test_aux = AUXSIGNAL.data;
    fe_test_aux = filtfilt(1/round(fs/10)*ones(round(fs/10),1), 1, fe_test_aux);
    
    % Load the full annotation to extract maximal excitation (only for
    % simulated data !)
    annotation_file = getfield(dir('*.ann'), 'name');
    fe_firings = ann2firings(annotation_file, fs);
    %fe_firings{numel(fe_firings)} = [];
    fe_spikes = firings2spikes(fe_firings, length(fe_test_emg));
        
    % Append to dataset
    test_emg = [test_emg; zeros(fe_min_pause * fs, size(fe_test_emg,2)); fe_test_emg];
    test_aux = [test_aux; zeros(fe_min_pause * fs, size(fe_test_aux,2)); fe_test_aux];
    
    test_spk = [test_spk, zeros(size(test_spk,1), size(fe_spikes,2) - size(test_spk,2))];
    test_spk = [test_spk; zeros(fe_min_pause * fs, size(test_spk,2))];
    test_spk = [test_spk; fe_spikes, zeros(size(fe_spikes,1), size(test_spk,2) - size(fe_spikes,2))];
   
    cd('..');
end

fprintf("\nTraining/Testing data loading complete! \n");
fsl = 64;
catch err
    fprintf("\nLoading data failed!");
    cd(current_path);
    rethrow(err);
end
    
cd(current_path);
clear fe_*