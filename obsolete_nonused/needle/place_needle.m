function [ needle_pts, needle_normals ] = place_needle(needle_pts, needle_normals, origin, orientation)

needle_pts = rodrigues_rot(needle_pts, [1, 0, 0],  orientation(1)); % roll
needle_pts = rodrigues_rot(needle_pts, [0, 1, 0],  orientation(2)); % pitch
needle_pts = rodrigues_rot(needle_pts, [0, 0, 1],  orientation(2)); % yaw

needle_normals = rodrigues_rot(needle_normals, [1, 0, 0],  orientation(1));
needle_normals = rodrigues_rot(needle_normals, [0, 1, 0],  orientation(2));
needle_normals = rodrigues_rot(needle_normals, [0, 0, 1],  orientation(3));

needle_pts = needle_pts + repmat(origin(:)', size(needle_pts,1), 1);

end

