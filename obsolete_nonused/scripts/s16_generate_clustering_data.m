
clusters_n = numel(detectable_ind); %Number of clusters to take
    
seg_half_len = floor(7.5/2/1000 * fs);
seg_full_len = 2*seg_half_len + 1;
seg_n = 100;
nmj_j = 25e-6;
max_std = 0;

%%
ap_centers_detectable = zeros(size(detectable_ind));
true_centers = zeros(clusters_n, seg_full_len);
for c = 1:clusters_n
    temp_ind = detectable_ind(c);
    
    % Get the non-jittered version of the muap
    clustering_data{c} = zeros(seg_n, seg_full_len);
    temp_muap = [MUs(temp_ind).calc_muap(0); zeros(max_muap_len-size(MUs(temp_ind).muap,1),size(MUs(temp_ind).muap,2))];
    temp_traj = traj_mixing_mat(0, traj_n_nodes, n_channels);   
    temp_muap = temp_muap * temp_traj' * electrode_diff_mat';
    
    
    % Get the max square of the energy (for additive noise scaling)
    if std(temp_muap) > max_std
        max_std = std(temp_muap);
    end
    
    % Generating the segments
    ap_centers_detectable(c) = find_ap_center_gauss(temp_muap);
    temp_muap = take_scope(temp_muap, ap_centers_detectable(c), seg_half_len); 
    true_centers(c, :) = temp_muap;
    
    for s = 1:seg_n
        temp_traj = traj_mixing_mat(s/seg_n, traj_n_nodes, n_channels);    
        temp_muap = [MUs(temp_ind).calc_muap(nmj_j); zeros(max_muap_len-size(MUs(temp_ind).muap,1),size(MUs(temp_ind).muap,2))];
        temp_muap = temp_muap * temp_traj' * electrode_diff_mat';
        temp_muap = take_scope(temp_muap, ap_centers_detectable(c), seg_half_len);
        
        clustering_data{c}(s,:) = temp_muap; 
    end
    fprintf('%d clusters data generated\n', c);
end

%% 
cd 'simulation_output';
mkdir 'clustering_data';
cd 'clustering_data';
save('clustering_data', 'clustering_data', 'clusters_n', 'seg_full_len', 'seg_n', 'nmj_j', 'max_std');

fid = fopen('commentary_on_clustering_data.txt', 'w');
% Preambule
fprintf(fid, 'These isolated EMG action potentials are simuluated using simulator created in Laboratoire des Sciences du Numerique de Nantes.\n\n');
fprintf(fid, 'Number of clusters: %d\n', clusters_n);
fprintf(fid, 'Neromuscular jitter standard deviation %d', nmj_j);
fprintf(fid, 'Segments per cluster %d', seg_n);
fprintf(fid, 'These segments contain jitter, but no additive noise, consider adding additive noise by, for example: segment = segment + max_std/100 * randn(size(segment)) ');
cd ..

clear seg_* C temp_*