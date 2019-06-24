%% Define MU parameters:
nmj_jitter = 35e-6;         % standard deviation, default value is 35e-6;
                        % Taken from Hamilton Stashuk

branch_v = [10000, 2000]; % NMJ branch conduction velocity; [mm/s] 
                         % 5000mm/s from Yann Pereon

                         
                         
%% MUs objects initialization
if ~exist('MUs', 'var')
    MUs(N) = MU_Sim;
end
for i = 1:N
    MUs(i) = MU_Sim(mf_centers(mf_assignment==i, :), 0, Lmuscle, ...
                    mf_diameters(mf_assignment == i), ...
                    mf_cv(mf_assignment==i), branch_v);
    
    MUs(i).nominal_center = mu_centers(i,:);
end



%% Generate neuromuscular junction coordinates distribution
endplate_area_center = Lmuscle/2;
n_branches = 1 + round(log(sz/sz(1)));

for i = 1:N
   arborization_z_std = 0.5 + sum(sz(1:i))/sum(sz) * 1.5;
   branches_z_std =     1.5 +  sum(sz(1:i))/sum(sz) * 4;
   MUs(i).sim_nmj_branches_two_layers(n_branches(i), endplate_area_center, branches_z_std, arborization_z_std);
    
%     disp('Attention, simplified nmj distribution is used')
%     branches_z_std = 1 + sum(sz(1:i))/sum(sz) * 2.5;
%     MUs(i).sim_nmj_branches_gaussian(endplate_area_center, branches_z_std);
end



%% Calculate MUs sfaps
for i = 1:N
    MUs(i).calc_sfaps(dt, dz, electrode_pts, electrode_normals);
    fprintf('%d sfaps generated\n', i);
end



%% Calculate MUAP templates (no jitter)
for i = 1:N
    MUs(i).calc_muap(0);
end



%% Graphics
%% Plot MUAPs
figure; ax =[]; max_muap_len = 0;
side = ceil(sqrt(N));
for i = 1:side
    for j = 1:side
        ind = (i-1)*side + j;   
        ax(ind) = subplot(side, side, ind);
        if ind <= N
            muap_to_plot = MUs(ind).muap;
            timeline = (0:size(muap_to_plot,1)-1)/fs * 1000; % ms
            plot(timeline, muap_to_plot * traj_mixing_mat(0, traj_n_nodes, n_channels)' * electrode_diff_mat');
        end
        axis tight;
        
        if numel(muap_to_plot) > max_muap_len
            max_muap_len = numel(muap_to_plot);
        end
    end
end
figure; 
for i = 1:side
    muap_to_plot = [];
    for j = 1:side
        ind = (i-1)*side + j; 
        if ind <= N
            muap_to_add = MUs(ind).muap * electrode_diff_mat';
            muap_to_plot = [muap_to_plot; muap_to_add(1:10/1000*fs,:)];
        end
    end
    subplot(side,1,i);
    timeline = (0:size(muap_to_plot,1)-1)/fs * 1000; % ms
    plot(timeline, muap_to_plot);
end
    
clear muap_to_plot muap_to_add timeline ax ind i j side
%linkaxes(ax(:));



%% Plot neuromuscular junctions z-coordinates distribtions
figure; ax =[]; set(gcf, 'Name', 'NMJ Z Coordinates Distribution, mm');
side = ceil(sqrt(N));
for i = 1:side
    for j = 1:side
        ind = (i-1)*side + j;
        ax(ind) = subplot(side, side, ind);
        if ind <= N
            histogram(MUs(ind).nmj_z)
            axis tight;
            xlim([min(0, min(MUs(ind).nmj_z)), max(Lmuscle, max(MUs(ind).nmj_z))]);
            set(gcf, 'position', [8 180 1433 530]);
        end
    end
end
clear ax side i j ind 


%% Plot axon branch lengths distribution
figure; ax =[]; set(gcf, 'Name', 'Branch Length Distribution, mm');
side = ceil(sqrt(N));
for i = 1:side
    for j = 1:side
        ind = (i-1)*side + j;
        ax(ind) = subplot(side, side, ind);
        if ind <= N, histogram(sum(MUs(ind).nerve_paths,2)), end
        axis tight;
    end
end
clear side ind i j

%% Plot axon propagation delays distribution
figure; ax =[]; set(gcf, 'Name', 'NMJ Delays Distribution, ms');
global_delays = [];
side = ceil(sqrt(N));
for i = 1:side
    for j = 1:side
        ind = (i-1)*side + j;
        ax(ind) = subplot(side, side, ind);
        if ind <= N
            histogram(MUs(ind).mnap_delays * 1000)
            axis tight;
            global_delays = [global_delays; MUs(ind).mnap_delays*1000];
        end
    end
end

figure; histogram(global_delays); 
xlabel('Delays between axon and it''s neuromuscular junctions, ms');
ylabel('Histogram of delays');

clear global_delays ind ax i j side


%% MUs' visibility by the electrode 
% ... in terms of muap energy per channel:
muaps_stds = zeros(N, n_channels);
for ch = 1:n_channels
    for mu = 1:N
        muaps_stds(mu, ch) = std(MUs(mu).muap * electrode_diff_mat(ch,:)');
    end
end
figure; 
subplot(1,2,1); plot(muaps_stds); 
title('MUs'' MUAPs energies seen by the electrode'); xlabel('MUs'); ylabel('STD');

% ... in terms of mean distance
fibers_distance = zeros(N, n_points);
for pt = 1:n_points
    for mu = 1:N
        dists = sqrt((MUs(mu).mf_centers(:,1) - electrode_pts(pt,1)).^2 + (MUs(mu).mf_centers(:,2) - electrode_pts(pt,2)).^2);
        fibers_distance(mu, pt) = mean(dists);
    end
end
subplot(1,2,2); plot(fibers_distance); 
title('Mean distances from MFs to the electrode');
xlabel('MUs'); ylabel('mm');

% ... in terms of mean electric distance (?)
% figure; hold on;
% for pt = 1:n_points
%     subplot(ceil(sqrt(n_points)), ceil(sqrt(n_points)), pt);
%     r = sqrt((mf_centers(:,1) - electrode_pts(pt,1)).^2 + (mf_centers(:,2) - electrode_pts(pt,2)).^2);
%     r(r<55e-3) = 55e-3;
%     ampls = get_elementary_current_response(0,0,r);
%     plot3(mf_centers(:,1), mf_centers(:,2), ampls, 'b.');
% end

% ... in terms of sfap power
% figure;
% for ch = 1:size(electrode_diff_mat,1)
%     subplot(ceil(sqrt(size(electrode_diff_mat,1))), ceil(sqrt(size(electrode_diff_mat,1))), ch); hold on;
%     for m = 1:N
%         sum_sfaps = zeros(size(squeeze(MUs(m).sfaps(:,1,:))));
%         for pt = 1:n_points
%             sum_sfaps = sum_sfaps + electrode_diff_mat(ch,pt) * squeeze(MUs(m).sfaps(:,pt,:));
%         end
%         plot3(MUs(m).mf_centers(:,1), MUs(m).mf_centers(:,2), std(sum_sfaps), 'b.');
%     end
% end


clear ampls r pt dists fibers_distance sum_sfaps ch



%% Jitter demonstration
figure; hold on;
n_to_plot = 22;
timeline = (0:numel(MUs(n_to_plot).calc_muap(0))-1)/fs*1000;
muap_to_plot = MUs(n_to_plot).calc_muap(0) * electrode_diff_mat';
plot(timeline(1:8/1000*fs), muap_to_plot(1:8/1000*fs, :) , 'k', 'linewidth', 2);
for i = 1:20
    muap_to_plot = MUs(n_to_plot).calc_muap(35e-6) * electrode_diff_mat';
    plot(timeline(1:8/1000*fs), muap_to_plot(1:8/1000*fs, :), 'g', 'linewidth', 0.2);
end
muap_to_plot = MUs(n_to_plot).calc_muap(0) * electrode_diff_mat';
plot(timeline(1:8/1000*fs), muap_to_plot(1:8/1000*fs, :) , 'k', 'linewidth', 2);

legend('Original MUAP', '20 APs with jitter');
%title('Jitter effect on the MUAP');
xlabel('Time, ms'); ylabel('Amplitude, arbitrary units');

clear  muap_to_plot timeline n_to_plot



%% Clear stuff
clear i j ind saved_scope x_circle y_circle muscle_border 
