%load('/Users/akmbpro/Nextcloud/EMG Data/simulated/physiological/Paper_On_Simulation/multichannel/multichannel_mdl.mat')

targetMU = 83; 
scope_time_ms = 10;
scope_time = 1:(20 * scope_time_ms);
scope_channels = 1:16;    

%% Analysis of how muap changes along electrode trajectory (requires single channel electrode);
upsample = 1;
snapshots = linspace(0,1, numel(electrode.traj_transforms(:,1) * upsample));

muap_z = zeros(size(MUs(targetMU).muap,1), length(snapshots), numel(scope_channels));
for c = 1 : length(snapshots)
    dummy = MUs(targetMU).muap * electrode.traj_mixing_mat(snapshots(c), electrode.n_nodes, electrode.n_channels)' * electrode.diff_mat';   
    for j = scope_channels
        muap_z(:,c,j) = dummy(:,j);
    end
end


%% Plot muscle and electrode
figure();
% Top view
l = line([0,0,Lmuscle, Lmuscle,0], [-Rmuscle, Rmuscle, Rmuscle, -Rmuscle, -Rmuscle], 'linewidth', 1); hold on;
text(6.9*Lmuscle/8, Rmuscle+0.5, 'Muscle border', 'fontsize', 8, 'color', 'k');
set(l, 'Color', 'k');
plot(electrode.pts_init(scope_channels, 3), electrode.pts_init(scope_channels, 1), 'k');
for c = scope_channels
    plot(electrode.pts_init(c, 3), electrode.pts_init(c, 1), '.', 'markersize', 10, 'color', [0.5, 0.5, 0.5]);
    %text(electrode.pts_init(c, 3)-0.5, electrode.pts_init(c, 1)+0.5, num2str(c), 'fontsize', 8);
    text(electrode.pts_init(c, 3), electrode.pts_init(c, 1)+0.5, num2str(c), 'fontsize', 8, 'horizontalalignment', 'center');
end

xlabel('Z, mm'); ylabel('Y, mm');
%axis equal
ylim([-6, 14]);
lborder = 28; rborder = 46;
xlim([lborder, rborder]);
smuap = max(abs(muap_z(:))); 
sx = 2; sy = 1.5;
%bx = linspace(lborder + 1, rborder - sx - 1, electrode.n_channels);
bx = electrode.pts_init(1:end-1, 3) + (electrode.pts_init(2, 3) - electrode.pts_init(1, 3))/2;
%by = Rmuscle + sy + ones(size(bx)); by(2:2:end) = 2*(sy+1) + by(2:2:end);
by0 = (Rmuscle + sy + 1) * ones(size(bx)); 
by = by0; by(2:3:end) = (sy+1) + by0(2:3:end); by(3:3:end) = 2*(sy+1) + by0(3:3:end);


for c = scope_channels(1:end-1)
    plot(bx(c) + sx*linspace(0,1,numel(scope_time)), by(c) + sy*muap_z(scope_time,1,c)/smuap, 'k'); %max(abs(muap_z(scope_time,1,c))));
    plot([bx(c), bx(c)], [by(c), by(c)+sy/2], 'k');
    plot(bx(c), by(c) + sy/2, 'ko', 'markersize', 1, 'MarkerFaceColor', 'k');
    text(bx(c) - 0.5, by(c) + sy/2 + 0.5, ['Ch ', num2str(c)], 'fontsize', 8, 'horizontalalignment', 'left');
    if c < 13
        cent = find_ap_center_gauss(muap_z(scope_time,1,c));
        %plot([bx(c), bx(c)]+sx*(cent/numel(scope_time)), by(c) + [0.8,1.2], '--', 'color', [0.2, 0.5, 0.3], 'linewidth', 1);
    end
    if ~mod(c+2,3) && (c < 15)
        %plot([electrode.pts_init(c, 3), electrode.pts_init(c, 3)], [electrode.pts_init(c, 1), by(c)], 'k:');
        plot([bx(c), bx(c)], [electrode.pts_init(c, 1) + (electrode.pts_init(2, 1)-electrode.pts_init(1, 1))/2, by(c)], 'k:');
    end
    
end

%% Plot innervation area
mfs = mu_pool.mf_centers(mu_pool.assignment == targetMU, :);
hull = convhull(mfs(:,1), mfs(:,2));
terr_top = max(mfs(:,1)-0.5);
terr_bot = min(mfs(:,1)+0.5);

patch([0, Lmuscle, Lmuscle, 0],  [terr_top,terr_top,terr_bot,terr_bot], 'b', 'facealpha', 0.25, 'edgealpha', 0);
%plot([0,Lmuscle], [terr_top,terr_top], '--b');
%plot([0,Lmuscle], [terr_bot,terr_bot], '--b');
text(6.9*Lmuscle/8, terr_top+0.5, 'MU territory', 'fontsize', 8, 'color', 'b');

%% Plot time annotation
plot([bx(1), bx(1)+sx], [by(1),by(1)]-1 , 'k'); %max(abs(muap_z(scope_time,1,c))));
plot([bx(1), bx(1)], [by(1)-0.2,by(1)+0.2]-1 , 'k'); %max(abs(muap_z(scope_time,1,c))));
plot([bx(1), bx(1)]+sx, [by(1)-0.2,by(1)+0.2]-1 , 'k'); %max(abs(muap_z(scope_time,1,c))));
text(bx(1) + sx/2, by(1)-1.75, [num2str(scope_time_ms), ' ms'], 'fontsize', 8, 'horizontalalignment', 'center');

%% Plot axes
%quiver(lborder+1,electrode.pts(4,1),1,0, 'linewidth', 1, 'maxheadsize', 1.5, 'Color', 'k');
%quiver(lborder+1,electrode.pts(4,1),0,1.75, 'linewidth', 1, 'maxheadsize', 0.5, 'Color', 'k');
ax = 0.15; ay = 0.4; axx = 0.05; ayy = 0.07;
annotation('arrow', [ax ax+axx],  [ay ay], 'color', 'k', 'linewidth', 0.5, 'headwidth', 5, 'headlength', 5);
annotation('arrow', [ax ax],  [ay ay+ayy], 'color', 'k', 'linewidth', 0.5, 'headwidth', 5, 'headlength', 5);
annotation('textbox', [ax+axx, ay, 0.025, 0.035], 'string', 'z', 'color', 'k', 'edgecolor', 'none', 'fontsize', 8);
annotation('textbox', [ax, ay+ayy, 0.025, 0.025], 'string', 'x', 'color', 'k', 'edgecolor', 'none', 'fontsize', 8);

%%
set(gca, 'ytick', []);
color = get(gcf,'Color');
set(gca,'YColor',color,'TickDir','out')

%% Render
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', [0, 0, 14.8, 10.5]); %A8 [0, 0, 7.4, 5.2], A7 [0, 0, 10.5, 7.4]); A6 [0, 0, 14.8, 10.5]
set(gcf, 'papersize', [21 29.7]); %A4
set(gca, 'fontunits', 'points', 'fontsize', 8);

print -dpng -r800 multichannel_muaps.png
%print -depsc2 multichannel_muaps.eps


