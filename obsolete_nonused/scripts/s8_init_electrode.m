%% Point electrode
% electrode_type = 'single point';
% electrode_pts_init = [0, 0, 2/3*Lmuscle];
% electrode_normals = zeros(size(electrode_pts_init));
% simplified_pts = electrode_pts_init;
% 
% electrode_diff_mat = 1;
% n_electrodes = size(electrode_diff_mat,2); 
% n_channels = size(electrode_diff_mat,1);
% n_points = size(electrode_pts_init,1);
% 

%% Or define a point differential electrode
electrode_type = '2-point differential'; electrode_type_short = '2p_diff';
electrode_step = 1; %mm
electrode_pts_init = [1, 0, 3/4*Lmuscle - electrode_step/2;
                      1, 0, 3/4*Lmuscle + electrode_step/2];
electrode_normals = zeros(size(electrode_pts_init));
simplified_pts = electrode_pts_init;

electrode_diff_mat = [-1, 1];
n_electrodes = size(electrode_diff_mat,2); 
n_channels = size(electrode_diff_mat,1);
n_points = size(electrode_pts_init,1);


%% Or define a needle
% electrode_type = 'needle'; electrode_type_short = 'needle';
% needle_len = 15; needle_out_rad = 0.3; needle_in_rad = 0.15;
% needle_roll = -pi; needle_pitch = 0; needle_yaw = 0;
% needle_x = 0; needle_y = 0; needle_z = 2/3*Lmuscle;
% 
% [electrode_pts_init, electrode_normals, electrode_signs, simplified_pts] = ...
%                 get_needle(needle_len, needle_out_rad, needle_in_rad, 6, needle_len, 5);
% 
% [electrode_pts_init, electrode_normals] = ...
%                 place_needle(electrode_pts_init, electrode_normals, [needle_x,needle_y,needle_z], [needle_roll, needle_pitch, needle_yaw]);
% 
% [simplified_pts] = ...
%                 place_needle(simplified_pts, [0,0,0], [needle_x,needle_y,needle_z], [needle_roll, needle_pitch, needle_yaw]);
% 
% % Filter out the points that are outside the muscle
% electrode_normals(sqrt(sum(electrode_pts_init(:, 1:2).^2, 2)) > Rmuscle, :) = [];
% electrode_signs(sqrt(sum(electrode_pts_init(:, 1:2).^2, 2)) > Rmuscle) = [];
% electrode_pts_init(sqrt(sum(electrode_pts_init(:, 1:2).^2, 2)) > Rmuscle, :) = [];
% 
% 
% electrode_diff_mat = electrode_signs(:)';
% n_electrodes = size(electrode_diff_mat,2); 
% n_channels = size(electrode_diff_mat,1);
% n_points = size(electrode_pts_init,1);

%% Or a multichannel iEMG
% electrode_type = 'intramuscular_array'; electrode_type_short = 'array';
% electrode_step = 1;
% needle_roll = -pi/4; needle_pitch = 0; needle_yaw = 0;
% %needle_roll = 0; needle_pitch = -pi/2; needle_yaw = 0;
% needle_x = 0; needle_y = 1; needle_z = 2/3*Lmuscle;
% %needle_x = 1; needle_y = - Rmuscle; needle_z = 2/3*Lmuscle;
% 
% %electrode_pts_init = [zeros(3,2), (0:1:2)'*electrode_step];
% electrode_pts_init = [zeros(5,2), (0:1:4)'*electrode_step];
% %electrode_pts_init = [zeros(9,2), (0:1:8)'*electrode_step];
% %electrode_pts_init = [zeros(17,2), (0:1:16)'*electrode_step];
% %electrode_pts_init = [zeros(129,2), (0:1:128)'*electrode_step];
% electrode_pts_init = rodrigues_rot(electrode_pts_init, [1, 0, 0],  needle_roll);
% electrode_pts_init = rodrigues_rot(electrode_pts_init, [0, 1, 0],  needle_pitch);
% electrode_pts_init = rodrigues_rot(electrode_pts_init, [0, 0, 1],  needle_yaw);
% electrode_pts_init = electrode_pts_init + repmat([needle_x, needle_y, needle_z], size(electrode_pts_init,1), 1);
% 
% electrode_normals = zeros(size(electrode_pts_init,1),3);
% electrode_signs = ones(size(electrode_pts_init,1),1);
% simplified_pts = electrode_pts_init;
% 
% electrode_diff_mat = eye(size(electrode_pts_init,1)-1,size(electrode_pts_init,1)) - circshift(eye(size(electrode_pts_init,1)-1,size(electrode_pts_init,1)), 1, 2);
% n_electrodes = size(electrode_diff_mat,2); 
% n_channels = size(electrode_diff_mat,1);
% n_points = size(electrode_pts_init,1);
% electrode_type = [electrode_type, num2str(n_channels),'ch'];
% electrode_type_short = [electrode_type_short, num2str(n_channels),'ch'];

%% Define electrode trajectory
% Simple case: two-point line
% Electrode is modeled as a rigid point cloud, so that its own geometry does
% not change while moving (which can be wrong in case of fine-wire
% electrodes, so pay attention to that).
% The trajectory is described as a number o transformations from
% the initial node to the current one, it contains six numbers: x, y, z
% and alpha, beta, gamma (yaw, pitch and roll angles).

% Stationary electrode:
traj_transforms =    [0,0,0,0,0,0];

%% Slight shifting along x axis (transversally to the fibers)
%traj_transforms = linspace(0,electrode_step,10);
%traj_transforms = [traj_transforms(:), zeros(length(traj_transforms), 5)];

%% Slight shifting along z axis (longitudinally to the fibers)
%traj_transforms = linspace(0, electrode_step, 10);
%traj_transforms = [zeros(length(traj_transforms), 2), traj_transforms(:), zeros(length(traj_transforms), 3)];

%% Scan across muscle
%traj_transforms = linspace(-Rmuscle, Rmuscle , round(2*Rmuscle/0.1));
%traj_transforms = [traj_transforms(:), zeros(length(traj_transforms), 5)];

% traj_transforms =    [0.00,0,0,0,0,0;
%                       0.10,0,0,0,0,0;
%                       0.20,0,0,0,0,0;
%                       0.30,0,0,0,0,0;
%                       0.40,0,0,0,0,0;
%                       0.50,0,0,0,0,0]/2; %x,y,z, yaw,pitch,roll;

% Other trajectories:

% All observation points along the trajectory:
traj_n_nodes = size(traj_transforms,1);
electrode_pts = [];
for i = 1:traj_n_nodes
    new_node = rotate_and_translate(electrode_pts_init, traj_transforms(i, 4:6), traj_transforms(i, 1:3));
    electrode_pts = [electrode_pts; new_node];
end

% Extended electrode matrix
electrode_diff_mat = repmat(electrode_diff_mat, 1, traj_n_nodes);

% Position-dependent mixing matrix to approximate the simulation output 
% between the nodes.
% Parameter t is bounded between 0 (initial position) and 1 (finish of the
% trajectory). You can map any behaviour into that parameter, for example, 
% the contraction force, or time. 

%% What makes the electrode translate: time or force.
% Time (normalized, relative to simulation duration)
traj_mixing_fun = @(t, n_nodes, node) max(0, 1 - (n_nodes-1) * abs(t - (node-1)/(n_nodes-1+eps)) );
traj_mixing_mat = @(t, n_nodes, n_channels) diag(reshape(repmat(traj_mixing_fun(t, n_nodes, 1:n_nodes), n_points, 1), [], 1));

% Force (normalized to maximum value)
%traj_parameter_map = force/max(force); % Force-dependent translation of the electrode
%traj_parameter_map_type = 'Force-dependent';


%% Show stuff
%% Plot muscle and territories
x_circle = linspace(-Rmuscle, Rmuscle, 100); y_circle = sqrt(Rmuscle^2 - x_circle.^2);
x_circle = [x_circle, x_circle(end:-1:1)]; y_circle = [y_circle, -y_circle(end:-1:1)];
figure();

% Cross-section view
subplot(2,2,1);
plot(x_circle, y_circle, 'k'); hold on;
for i = 1:N
    text(mu_centers(i,1), mu_centers(i,2), num2str(i));
end
axis equal
axis manual
saved_scope = axis;
plot(simplified_pts(:,1), simplified_pts(:,2), 'b*');
for i = 1:n_electrodes
    plot(electrode_pts(i:n_electrodes:end, 1), electrode_pts(i:n_electrodes:end, 2), 'b');
end
title('Cross-section view'); xlabel('Height, mm'); ylabel('Width, mm');

% Top view
subplot(2,2,3); hold on;
l = line([-Rmuscle, Rmuscle, Rmuscle, -Rmuscle, -Rmuscle], [0,0,Lmuscle, Lmuscle,0]);
set(l, 'Color', 'k');
for i = 1:N
    line([mu_centers(i,1), mu_centers(i,1)], [0, Lmuscle]);
end

plot(simplified_pts(:,1), simplified_pts(:,3), '*');
for i = 1:n_electrodes
    plot(electrode_pts(i:n_electrodes:end, 1), electrode_pts(i:n_electrodes:end, 3), 'b');
end
xlim(saved_scope([1,2]));
title('Top view'); xlabel('Width, mm'); ylabel('Length, mm');

% Left view
subplot(2,2,2); hold on;
l = line([0, Lmuscle, Lmuscle, 0, 0], [-Rmuscle, -Rmuscle, Rmuscle, Rmuscle, -Rmuscle]);
set(l, 'Color', 'k');
for i = 1:N
    line([0, Lmuscle], [mu_centers(i,2), mu_centers(i,2)]);
end
plot(simplified_pts(:,3), simplified_pts(:,2), '*');
for i = 1:n_electrodes
    plot(electrode_pts(i:n_electrodes:end, 3), electrode_pts(i:n_electrodes:end, 2), 'b');
end
title('Left view'); xlabel('Length, mm'); ylabel('Width, mm');

% MF innervation
% subplot(2,2,4);
% plot(x_circle, y_circle, 'k'); hold on;
% for i = N:-1:1
%     plot(MUs(i).mf_centers(:,1), MUs(i).mf_centers(:, 2), '.', 'MarkerSize', 5);
% end
% axis equal
% title('MU assignment');

clear x_circle y_circle l

