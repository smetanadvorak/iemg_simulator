analysisN = 17;
to_take_time = 1:175;
to_take_traj = 1:35;

%% Analysis of how muap changes along electrode trajectory (requires single channel electrode);
if electrode.n_channels > 1
    error('Number of channels shouldn''t be larger than 1');
end

%% Draw MUAPs in the nodes
snapshots = 0 : 1/numel(electrode.traj_transforms(:,1)) : 1;
muap_z_skeleton = zeros(size(MUs(analysisN).muap,1), length(snapshots));

for i = 1 : length(snapshots)
    muap_z_skeleton(:,i) = MUs(analysisN).muap * electrode.traj_mixing_mat(snapshots(i), electrode.n_nodes, electrode.n_channels)' * electrode.diff_mat';   
end

figure; hold on;
muap_z_skeleton = muap_z_skeleton(to_take_time,:);
[X,Y] = meshgrid(snapshots*(electrode.n_nodes-1)*electrode.traj_step, (1:size(muap_z_skeleton,1))/fs * 1000);

for i = 1:electrode.n_nodes
    plot3(X(:,i), Y(:,i), muap_z_skeleton(:,i), 'k', 'linewidth',1);
end

%% Draw muaps in between the nodes (interpolated)
upsample = 2;
snapshots = 0 : 1/numel(electrode.traj_transforms(:,1))/upsample : 1;
to_take_traj = to_take_traj(1) : to_take_traj(end) * upsample;

muap_z = zeros(size(MUs(analysisN).muap,1), length(snapshots));

for i = 1 : length(snapshots)
    muap_z(:,i) = MUs(analysisN).muap * electrode.traj_mixing_mat(snapshots(i), electrode.n_nodes, electrode.n_channels)' * electrode.diff_mat';   
end

muap_z = muap_z(to_take_time,:);
[X,Y] = meshgrid(snapshots()*(electrode.n_nodes-1)*electrode.traj_step, (1:size(muap_z,1))/fs * 1000);
%[X,Y] = meshgrid(electrode.traj_transforms(:,1), (1:size(muap_z,1))/fs * 1000);

surf(X,Y,muap_z);
alpha 0.50

ylabel('Time, ms'); 
%xlabel('Trajectory nodes');
xlabel('Translation along trajectory, mm');
zlabel('Amplitude, arbitrary units');
%figure2page;
%export_fig MUAP_trajectory.png -r300;
%saveas(gcf, 'MUAP_trajectory', 'png');

%% Draw noise level
[X,Y] = meshgrid(snapshots()*(electrode.n_nodes-1)*electrode.traj_step, (1:size(muap_z,1))/fs * 1000);
%surf(X,Y, 3 * emg_noise_std*ones(size(X)));
alpha 0.50

