

len = 30;
out_rad = 1;
in_rad = 0.5;
circ_segments = 6;
lin_segments = 10;
tip_segments = 5;


[pts, normals, signs] = get_needle(len, out_rad, in_rad, circ_segments, lin_segments, tip_segments);

figure; 
plot3(pts(:,1), pts(:,2), pts(:,3), '*'); hold on;
plot3(pts(1,1), pts(1,2), pts(1,3), 'r*'); hold on;
plot3(pts(2,1), pts(2,2), pts(2,3), 'g*'); hold on;
quiver3(pts(:,1), pts(:,2), pts(:,3), normals(:,1), normals(:,2), normals(:,3), 0);

[pts, normals] = place_needle(pts, normals, [0,0,1], 0, pi/6, 0);

plot3(pts(:,1), pts(:,2), pts(:,3), 'k*'); hold on;
quiver3(pts(:,1), pts(:,2), pts(:,3), normals(:,1), normals(:,2), normals(:,3), 0, 'k');

axis equal