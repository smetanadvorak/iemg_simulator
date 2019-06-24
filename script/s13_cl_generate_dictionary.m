
%% Generate full dictionary
% MUAPs matrices are (samples x channel x motor_unit)
max_muap_len = 0;
for m = 1:mu_pool.N
    if size(MUs(m).muap,1) > max_muap_len
        max_muap_len = size(MUs(m).muap,1);
    end
end


% Full dictionary of forms in trajectory nodes
dictionary_full_traj = cell(electrode.n_nodes, 1);
for k = 1:electrode.n_nodes
    dictionary_full_traj{k} = zeros(max_muap_len, electrode.n_channels, mu_pool.N);
    for m = 1:mu_pool.N
        temp_muap = [MUs(m).muap; zeros(max_muap_len-size(MUs(m).muap,1),size(MUs(m).muap,2))];
        temp_traj = electrode.traj_mixing_mat((k-1)/(electrode.n_nodes-1+eps), electrode.n_nodes, electrode.n_channels);
        dictionary_full_traj{k}(:,:,m) = temp_muap * temp_traj' * electrode.diff_mat';
    end
end


% Full dictionary of initial forms (MUAPs in the intitial electrode position)
dictionary_full_init = dictionary_full_traj{1};

% Generate indices of MUAPs actually present in the signal (active motor
% units)
largest_active_full_mu = find(any(spikes), 1, 'last');


%% Generate detectable dictionary
% Full dictionary of forms in trajectory nodes
dictionary_detectable_traj = cell(electrode.n_nodes, 1);
for k = 1:electrode.n_nodes
    dictionary_detectable_traj{k} = zeros(max_muap_len, electrode.n_channels, numel(detectable_ind));
    for m = 1:numel(detectable_ind)
        temp_ind = detectable_ind(m);
        temp_muap = [MUs(temp_ind).muap; zeros(max_muap_len-size(MUs(temp_ind).muap,1),size(MUs(temp_ind).muap,2))];
        temp_traj = electrode.traj_mixing_mat((k-1)/(electrode.n_nodes-1+eps), electrode.n_nodes, electrode.n_channels);
        dictionary_detectable_traj{k}(:,:,m) = temp_muap * temp_traj' * electrode.diff_mat';
    end
end


% Detectable dictionary of initial forms (MUAPs in the intitial electrode position)
dictionary_detectable_init = dictionary_detectable_traj{1};


largest_active_detectable_mu = find(detectable_ind <= find(any(spikes),1,'last'), 1, 'last');

comment_full_init = 'dictionary_full_init(t, channel, motor unit)';
comment_full_traj = 'dictionary_full_traj{trajectory node}(t, channel, motor unit)';
comment_detectable_init = 'dictionary_detectable_init(t, channel, motor unit)';
comment_detectable_traj = 'dictionary_detectable_init{trajectory node}(t, channel, motor unit)';



%% In case of single channel electrode, automatically squeeze the dictionaries:
dictionary_detectable_init = squeeze(dictionary_detectable_init);
dictionary_full_init = squeeze(dictionary_full_init);

for i = 1:length(dictionary_detectable_traj)
    dictionary_detectable_traj{i} = squeeze(dictionary_detectable_traj{i});
end

for i = 1:length(dictionary_full_traj)
    dictionary_full_traj{i} = squeeze(dictionary_full_traj{i});
end

%%

% dictionary_full_init = zeros(max_muap_len, electrode.n_channels, mu_pool.N);
% for m = 1:mu_pool.N
%     temp_muap = [MUs(m).muap; zeros(max_muap_len-size(MUs(m).muap,1),size(MUs(m).muap,2))];
%     temp_traj =  electrode.traj_mixing_mat(0, electrode.n_nodes, electrode.n_channels);
%     dictionary_full_init(:,:,m) = temp_muap * temp_traj' * electrode.diff_mat';
% end

% dictionary_detectable_init = zeros(max_muap_len, electrode.n_channels, numel(detectable_ind));
% for m = 1:numel(detectable_ind)
%     temp_ind = detectable_ind(m);
%     temp_muap = [MUs(temp_ind).muap; zeros(max_muap_len-size(MUs(temp_ind).muap,1),size(MUs(temp_ind).muap,2))];
%     temp_traj = electrode.traj_mixing_mat(0, electrode.n_nodes, electrode.n_channels);
%     dictionary_detectable_init(:,:,m) = temp_muap * electrode.diff_mat';
% end

