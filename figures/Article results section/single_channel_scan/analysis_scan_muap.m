%load('/Users/akmbpro/Nextcloud/EMG Data/simulated/physiological/Paper_On_Simulation/single_channel_scan/simulation_mdl.mat')

targetMU = 77; %76
scope_time = 1:(20 * 8);
scope_traj = 1:35;

%% Analysis of how muap changes along electrode trajectory (requires single channel electrode);
if electrode.n_channels > 1
    error('Number of channels shouldn''t be larger than 1');
end

%% Draw MUAPs in the nodes
snapshots = 0 : 1/numel(electrode.traj_transforms(:,1)) : 1;

muap_z_skeleton = zeros(size(MUs(targetMU).muap,1), length(snapshots));
for i = 1 : length(snapshots)
    muap_z_skeleton(:,i) = MUs(targetMU).muap * electrode.traj_mixing_mat(snapshots(i), electrode.n_nodes, electrode.n_channels)' * electrode.diff_mat';   
end

figure; hold on;
muap_z_skeleton = muap_z_skeleton(scope_time,:);
[Xs,Ys] = meshgrid(snapshots*(electrode.n_nodes-1)*electrode.traj_step, (1:size(muap_z_skeleton,1))/fs * 1000);

for i = 1:electrode.n_nodes
    plot3(Xs(:,i), Ys(:,i), muap_z_skeleton(:,i), 'k', 'linewidth',0.5);
end

%% Draw muaps in between the nodes (interpolated)
upsample = 10;
snapshots = 0 : 1/numel(electrode.traj_transforms(:,1))/upsample : 1;
snapshots = snapshots(1:end-upsample);
scope_traj = scope_traj(1) : scope_traj(end) * upsample;
muap_z = zeros(size(MUs(targetMU).muap,1), length(snapshots));

for i = 1 : length(snapshots)
    muap_z(:,i) = MUs(targetMU).muap * electrode.traj_mixing_mat(snapshots(i), electrode.n_nodes, electrode.n_channels)' * electrode.diff_mat';   
end
muap_z = muap_z(scope_time,:);
[X,Y] = meshgrid(snapshots()*(electrode.n_nodes-1)*electrode.traj_step, (1:size(muap_z,1))/fs * 1000);

surf(X,Y,muap_z, 'edgecolor','none');
alpha 0.50
view([51.6000   32.4000]);
ylabel('Time, ms'); 
xlabel('Translation along trajectory, mm');
zlabel('Amplitude, arbitrary units');
grid on

%% Draw trajectory across territory
% Need: muscle cross section, territory polygone, scaling matrix, mfs
% positions
mfs = mu_pool.mf_centers(mu_pool.assignment == targetMU, :);
hull = convhull(mfs(:,1), mfs(:,2));
border = mu_pool.muscle_border;

% Normalize
mfs = mfs/max(border(:))/2;
border = border/max(border(:))/2;

% Fit to current coordinates
Z = zlim;
border(:,1) = border(:,1) * abs(min(Xs(:)) - max(Xs(:)));
border(:,2) = border(:,2) * abs(min(Z(:)) - max(Z(:)));
border(:,1) = border(:,1) + (min(Xs(:)) + max(Xs(:)))/2;
border(:,2) = border(:,2) + (min(Z(:)) + max(Z(:)))/2;

mfs(:,1) = mfs(:,1) * abs(min(Xs(:)) - max(Xs(:)));
mfs(:,2) = mfs(:,2) * abs(min(Z(:)) - max(Z(:)));
mfs(:,1) = mfs(:,1) + (min(Xs (:)) + max(Xs(:)))/2;
mfs(:,2) = mfs(:,2) + (min(Z(:)) + max(Z(:)))/2;

plot3(border(:,1), repmat(max(Y(:)), size(border,1), 1), border(:,2), 'k--', 'linewidth',1);
plot3(mfs(hull,1), repmat(max(Y(:)), numel(hull), 1), mfs(hull,2), 'b--', 'linewidth',0.5);
plot3(mfs(:,1), repmat(max(Y(:)), numel(mfs(:,1)), 1), mfs(:,2), '.b', 'markersize', 1);
plot3(Xs(1,:), max(Ys), zeros(size(Xs,2),1), 'ko', 'markerfacecolor', 'k', 'markersize', 1.5);

[~, max_ind] = max(mfs(hull, 2));
text(mfs(hull(max_ind),1), max(Y(:)), mfs(hull(max_ind),2)+1.5, 'MU territory', 'color', 'b', 'fontsize', 8, 'linewidth',1, 'horizontalalignment', 'right');
[~, max_ind] = max(border(:, 2));
text(border(max_ind,1), max(Y(:)), border(max_ind,2)+1.5, 'Muscle border', 'fontsize', 8, 'linewidth',1, 'horizontalalignment', 'right');



set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', [0, 0, 14.8, 10.5]); %A8 [0, 0, 7.4, 5.2], A7 [0, 0, 10.5, 7.4]); A6 [0, 0, 14.8, 10.5]
set(gcf, 'papersize', [21 29.7]); %A4
set(gca, 'fontunits', 'points', 'fontsize', 8);

%align_axislabel(1,gca);
% h = rotate3d;
% set(h,'ActionPreCallback',...
%     'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
% set(h,'ActionPostCallback',...
%     'set(gcf,''windowbuttonmotionfcn'','''')')

print -dpng -r800 singlechannel_scan_muap.png



% %% Draw noise level
% [X,Y] = meshgrid(snapshots()*(electrode.n_nodes-1)*electrode.traj_step, (1:size(muap_z,1))/fs * 1000);
% %surf(X,Y, 3 * emg_noise_std*ones(size(X)));
% alpha 0.50

