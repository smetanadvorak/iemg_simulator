N_analysis = 33;
to_take = 20:150;

%% Analysis of how muap changes along electrode trajectory (requires single channel electrode);
upsample = 0.5;
snapshots = linspace(0,1, numel(electrode.traj_transforms(:,1) * upsample));

muap_z = zeros(size(MUs(N_analysis).muap,1), length(snapshots), electrode.n_channels);
for i = 1 : length(snapshots)
    dummy = MUs(N_analysis).muap * electrode.traj_mixing_mat(snapshots(i), electrode.n_nodes, electrode.n_channels)' * electrode.diff_mat';   
    for j = 1:electrode.n_channels
        muap_z(:,i,j) = dummy(:,j);
    end
end

muap_z = muap_z(to_take,:,:);

% Channel 1
f1 = figure; subplot(1,2,1);
[X,Y] = meshgrid(snapshots, (1:size(muap_z(:,:,1),1))/fs * 1000);
%[Xp,Yp] = meshgrid(snapshots*(traj_n_nodes-1), (1:size(muap_z(:,:,1),1))/fs * 1000);
%[X,Y] = meshgrid(traj_transforms(:,1), (1:size(muap_z,1))/fs * 1000);
surf(X,Y,muap_z(:,:,1)); hold on;
alpha 0.1
ylabel('$t$, ms'); 
%xlabel('Trajectory nodes');
xlabel('$\lambda$');
zlabel('Amplitude, arbitrary units');
%align_axislabel([], gca)
 
% Channel 2
%f2 = figure; 
subplot(1,2,2);
[X,Y] = meshgrid(snapshots, (1:size(muap_z(:,:,2),1))/fs * 1000);
%[X,Y] = meshgrid(traj_transforms(:,1), (1:size(muap_z,1))/fs * 1000);
surf(X,Y,muap_z(:,:,2), 'linewidth', 0.25); hold on;
alpha 0.1
ylabel('$t$, ms'); 
%xlabel('Trajectory nodes');
xlabel('$\lambda$');
zlabel('Amplitude, arbitrary units');
%align_axislabel([], gca)


%% Draw MUAPs in the nodes
%upsample = 0.5;
%snapshots = linspace(0,1, numel(traj_transforms(:,1))*upsample);

% muap_z_skeleton = zeros(size(MUs(N_analysis).muap,1), length(snapshots), n_channels);
% for i = 1 : length(snapshots)
%     dummy = MUs(N_analysis).muap * traj_mixing_mat(snapshots(i), traj_n_nodes, n_channels)' * electrode_diff_mat';   
%     for j = 1:n_channels
%         muap_z_skeleton(:,i,j) = dummy(:,j);
%     end
% end
% muap_z_skeleton = muap_z_skeleton(to_take,:,:);

%Channel 1
%figure(f1);
subplot(1,2,1);
% [X,Y] = meshgrid(snapshots*(traj_n_nodes-1), (1:size(muap_z_skeleton(:,:,1),1))/fs * 1000);
% for i = 1:traj_n_nodes
%     plot3(X(:,i), Y(:,i), muap_z_skeleton(:,i,1), 'k', 'linewidth',0.5);
% end
plot3(X(:,1), Y(:,1), muap_z_skeleton(:,1,1), 'r', 'linewidth',3);
ax1 = gca;
%Channel 2 
%figure(f2); 
subplot(1,2,2);
% [X,Y] = meshgrid(snapshots*(traj_n_nodes-1), (1:size(muap_z_skeleton(:,:,2),1))/fs * 1000);
% for i = 1:traj_n_nodes
%     plot3(X(:,i), Y(:,i), muap_z_skeleton(:,i,2), 'k', 'linewidth',0.5);
% end
plot3(X(:,end), Y(:,end), muap_z_skeleton(:,end,2), 'r', 'linewidth',3);
ax2 = gca;

linkprop([ax2, ax1], {'xlim', 'ylim', 'zlim', 'cameraposition'});
%mesh(X,Y,muap_z_skeleton, 'Linewidth', 1);
% ylim([60,160]);
% xlim([0,20]);

%saveas(gcf, 'MUAP_trajectory', 'png');

