function pt = rotate_and_translate( pt, rpy, d )
% Rotates a cloud point pt (Nx3 matrix) around origin using rpy (3x1
% vector) of roll, pitch and yaw angles, then translates the cloud to d (3x1 vector);

pt = rodrigues_rot(pt, [1, 0, 0],  rpy(1)); % roll
pt = rodrigues_rot(pt, [0, 1, 0],  rpy(2)); % pitch
pt = rodrigues_rot(pt, [0, 0, 1],  rpy(3)); % yaw

pt = pt + repmat(d(:)', size(pt,1), 1); % translation

end

