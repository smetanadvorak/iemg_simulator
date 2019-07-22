
%% Set the electrode

% End-to-end array
% electrode = IntramuscularArray(11, 1, 'consecutive');
% electrode.set_position([-Rmuscle,0,2*Lmuscle/3], [-pi/2, 0, -pi/2]); 
% electrode.set_linear_trajectory(0, 1); %Static

% One-step scanning array
%electrode = IntramuscularArray(16, 1, 'consecutive');
%electrode.set_position([-Rmuscle,0,Lmuscle/2+5], [-pi/6, 0, -pi/2]); 
%electrode.set_linear_trajectory(1, 1); %One-step pull out

% Single channel differential electrode
electrode = IntramuscularArray(2, 1);
electrode.set_position([-2.5,0,2*Lmuscle/3], [-pi/2,0,-pi/2]);

% Scanning electrode
%electrode = IntramuscularArray(2, 1);
%electrode.set_position([-Rmuscle,0,2*Lmuscle/3], [-pi/2, 0, -pi/2]); 
%electrode.set_linear_trajectory(2*Rmuscle-1, 19);

%% Plot muscle and territories
x_circle = linspace(-Rmuscle, Rmuscle, 1000); y_circle = sqrt(Rmuscle^2 - x_circle.^2);
x_circle = [x_circle, x_circle(end:-1:1)]; y_circle = [y_circle, -y_circle(end:-1:1)];
figure(); set(gcf, 'position', [192  117 1260 761]);

% Cross-section view
subplot(2,2,1);
plot(x_circle, y_circle, 'k'); hold on;
for i = 1:mu_pool.N
    text(mu_pool.mn_pool.centers(i,1), mu_pool.mn_pool.centers(i,2), num2str(i));
end
axis equal
axis manual
saved_scope = axis;

for i = 1:size(electrode.pts,1)
    plot(electrode.pts(i, 1), electrode.pts(i, 2), 'k.');
end
for i = 1:size(electrode.pts_init,1)
    plot(electrode.pts_init(i, 1), electrode.pts_init(i, 2), 'bo');
end
title('Cross-section view'); xlabel('Height, mm'); ylabel('Width, mm');

% Top view
subplot(2,2,3); hold on;
l = line([-Rmuscle, Rmuscle, Rmuscle, -Rmuscle, -Rmuscle], [0,0,Lmuscle, Lmuscle,0]);
set(l, 'Color', 'k');
% for i = 1:mu_pool.N
%     line([mu_pool.mn_pool.centers(i,1), mu_pool.mn_pool.centers(i,1)], [0, Lmuscle]);
% end

for i = 1:size(electrode.pts,1)
    plot(electrode.pts(i, 1), electrode.pts(i, 3), 'k.');
end
for i = 1:size(electrode.pts_init,1)
    plot(electrode.pts_init(i, 1), electrode.pts_init(i, 3), 'bo');
end
xlim(saved_scope([1,2]));
title('Top view'); xlabel('Width, mm'); ylabel('Length, mm');
axis equal

% Left view
subplot(2,2,2); hold on;
l = line([0, Lmuscle, Lmuscle, 0, 0], [-Rmuscle, -Rmuscle, Rmuscle, Rmuscle, -Rmuscle]);
set(l, 'Color', 'k');
% for i = 1:mu_pool.N
%     line([0, Lmuscle], [mu_pool.mn_pool.centers(i,2), mu_pool.mn_pool.centers(i,2)]);
% end

for i = 1:electrode.n_points
    plot(electrode.pts(i:electrode.n_points:end, 3), electrode.pts(i:electrode.n_points:end, 2), '*');
end
title('Left view'); xlabel('Length, mm'); ylabel('Width, mm');
axis equal
clear x_circle y_circle l

