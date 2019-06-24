%% To Do:
% Make a function out of it!
% To solve: access to number of motor neurons
% Add: profile loading if generate a new dataset

%%
% This script loads the dataset for intent estimator training and
% evaluation. 

% 1) Set the dataset_path, it should contain one file:
% virtual_subject_mdl.mat; and two folders: train_data and test_data

% 2) Test/train data folders each contain a number of folders with
% contraction trials;

% 3) Each contraction trial should start and finish  with at least one 
% second of no activity (rest), in order to ensure that a trial doesn't
% affect the next trial in any way. 

% 4) Each trial folder contains at least these files:
% simulation_signals.mat with EMGSIGNAL and AUXSIGNAL structures
% sim_annotation_centered_detectable.ann with the decomposition
% and sim_annotation_centered_full.ann with full annotation (used only for
% visualization).

%dataset_path = '/Users/akmbpro/Nextcloud/EMG Data/simulated/physiological/data_for_estimator_tests';
dataset_path = '/Users/akmbpro/Nextcloud/EMG Data/simulated/physiological/virtual_subject_test2';

vs_mdl_file = 'simulation_mdl.mat';
signals_file = 'simulation_signals.mat';
annotation_file = 'sim_annotation_centered_full.ann';
training_path = 'training_data';
testing_path = 'testing_data';

%cd(fileparts(mfilename('fullpath')));
current_path = pwd;

try

cd(dataset_path);
load simulation_mdl;
mf_mdl.init_quasistatic_e2f_f2e_models(mu_pool);
%% Load the virtual subject model (optional)
%load(vs_mdl_file);

%% Load training data
% Make list of trials in this folder
cd(training_path);
dir_list = dir;
dir_list = dir_list([dir_list.isdir] == 1 & ~startsWith({dir_list.name}, '.'));
dir_list = {dir_list.name};

train_emg = [];
train_aux = [];
train_spikes = [];
train_max_aux = [];
fe_min_pause = 1; %In seconds;

fprintf("\nLoading training data ...");
for d = 1:numel(dir_list)
    cd(dir_list{d});
    
    % Load the signal;
    load(signals_file);
    fs = EMGSIGNAL.rate;
    fe_train_emg = EMGSIGNAL.data;
    fe_train_aux = AUXSIGNAL.data;
    fe_train_aux = filtfilt(1/round(fs/10)*ones(round(fs/10),1), 1, fe_train_aux);
    
    % Load the full annotation to extract maximal excitation (only for
    % simulated data !)
    fe_firings = ann2firings(annotation_file, fs);
    fe_firings{mu_pool.N+1} = []; fe_firings(mu_pool.N+1) = []; % Go up to the number of motor neurons
    fe_spikes = firings2spikes(fe_firings, length(fe_train_emg));
        
    % Append to dataset
    train_emg = [train_emg; zeros(fe_min_pause * fs, size(fe_train_emg,2)); fe_train_emg];
    train_aux = [train_aux; zeros(fe_min_pause * fs, size(fe_train_aux,2)); fe_train_aux];
    train_spikes = [train_spikes; zeros(fe_min_pause * fs, size(fe_spikes,2)); fe_spikes];
   
    cd ..;
end

train_spikes = train_spikes(:, detectable_ind);
max_active_mu = find(any(fe_spikes), 1, 'last');
max_excitation = ( mu_pool.mn_pool.rt(max_active_mu) + mu_pool.mn_pool.rt(min(mu_pool.N, max_active_mu+1)) )/2;
active_ind = detectable_ind(mu_pool.mn_pool.rt(detectable_ind) < max_excitation);
MEC = max(train_aux); % Maximal estimable contraction;

cd ..

%% Load the testing data
cd(testing_path);
dir_list = dir;
dir_list = dir_list([dir_list.isdir] == 1 & ~startsWith({dir_list.name}, '.'));
dir_list = {dir_list.name};

test_emg = [];
test_aux = [];
test_spikes = [];
test_max_aux = [];
fe_min_pause = 1; %In seconds;

fprintf("\nLoading testing data ...");
for d = 1:numel(dir_list)
    cd(dir_list{d});
    
    % Load the signal;
    load(signals_file);
    fs = EMGSIGNAL.rate;
    fe_test_emg = EMGSIGNAL.data;
    fe_test_aux = AUXSIGNAL.data;
    fe_test_aux = filtfilt(1/round(fs/10)*ones(round(fs/10),1), 1, fe_test_aux);
    
    % Load the full annotation to extract maximal excitation (only for
    % simulated data !)
    fe_firings = ann2firings(annotation_file, fs);
    fe_firings{mu_pool.N+1} = []; fe_firings(mu_pool.N+1) = []; % Go up to the number of motor neurons
    fe_spikes = firings2spikes(fe_firings, length(fe_test_emg));
        
    % Append to dataset
    test_emg = [test_emg; zeros(fe_min_pause * fs, size(fe_test_emg,2)); fe_test_emg];
    test_aux = [test_aux; zeros(fe_min_pause * fs, size(fe_test_aux,2)); fe_test_aux];
    test_spikes = [test_spikes; zeros(fe_min_pause * fs, size(fe_spikes,2)); fe_spikes];
   
    cd('..');
end

test_spikes = test_spikes(:, detectable_ind);
fprintf("\nTraining/Testing data loading complete! \n");

catch err
    fprintf("\nLoading data failed!");
    cd(current_path);
    rethrow(err);
end
    
cd(current_path);

%% Downsample the datasets
fsl = 100;
%% Define downsamplers
downsampler_aux = Downsampler(fs/fsl);
downsampler_spk = Downsampler(fs/fsl, size(train_spikes,2));
downsampler_act = Downsampler(fs/fsl, size(train_spikes,2));

%% Downsample the training data using online downsampler
train_aux_dwns = zeros(size(train_aux));
train_act_dwns = zeros(size(train_spikes));
train_spk_dwns  = zeros(size(train_spikes));
train_act = spikes2activity(train_spikes, length(fe_window)/fsl*fs);

for t = 1:size(train_spikes,1)
    train_aux_dwns(t) = downsampler_aux.update(train_aux(t));
    train_act_dwns(t,:) = downsampler_act.update(train_act(t,:));
    train_spk_dwns(t,:) = downsampler_spk.update(train_spikes(t,:));
end

train_aux_dwns(isnan(train_aux_dwns(:,1)),:) = [];
train_act_dwns(isnan(train_act_dwns(:,1)),:) = [];
train_spk_dwns(isnan(train_spk_dwns(:,1)),:) = [];

train_act_dwns = double(train_act_dwns >= 0.5);
train_spk_dwns = double(train_spk_dwns > 0);

%% Downsample the testing data
downsampler_aux.reset();
downsampler_spk.reset();
downsampler_act.reset();

test_aux_dwns = zeros(size(test_aux));
test_act_dwns = zeros(size(test_spikes));
test_spk_dwns  = zeros(size(test_spikes));
test_act = spikes2activity(test_spikes, length(fe_window)/fsl*fs);

for t = 1:size(test_spikes,1)
    test_aux_dwns(t) = downsampler_aux.update(test_aux(t));
    test_act_dwns(t,:) = downsampler_act.update(test_act(t,:));
    test_spk_dwns(t,:) = downsampler_spk.update(test_spikes(t,:));
end

test_aux_dwns(isnan(test_aux_dwns(:,1)),:) = [];
test_act_dwns(isnan(test_act_dwns(:,1)),:) = [];
test_spk_dwns(isnan(test_spk_dwns(:,1)),:) = [];

test_act_dwns = double(test_act_dwns >= 0.5);
test_spk_dwns = double(test_spk_dwns > 0);
%%
clear signals_file annotation_file fe_* EMGSIGNAL AUXSIGNAL current_path dir_list training_path testing_path