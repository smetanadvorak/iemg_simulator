%% Muscle geometry
% Everything is in millimeters. Muscle in vicinity of the electrode is approximated as a cylinder with
Lmuscle = 40; % [mm]
Dmf = 400; % Density of muscle fibers per square millimetre (Hamilton-Wright 2005)
Nmf = round((Rmuscle^2) * pi * Dmf); % Expected number of muscle fibers in the muscle. 


%% Muscle fibers coordinates generation
cd ./voronoi_sampling
mf_centers = csvread('./voronoi_pi1e5.csv');
cd ..   

mf_centers = (mf_centers - 5)/4;
[dists, inds] = sort(sqrt(mf_centers(:,1).^2 + mf_centers(:,2).^2));
mf_centers = mf_centers(inds(1:Nmf), :) / dists(Nmf+1)*Rmuscle;



%% MUs innervation numbers
innerv_num = round(innerv_area * Dmf);


%% Plot fibers distribution in muscle cross-section
phi_circle = linspace(0, 2*pi, 10000)'; %Should be sampled as finely as possible;
phi_circle = phi_circle(1:end-1); % No overlapping
muscle_border = [Rmuscle * cos(phi_circle), Rmuscle * sin(phi_circle)];
figure; plot(muscle_border, 'k');
plot(mf_centers(:,1), mf_centers(:,2), 'k.');

clear phi_circle




%% Assigning muscle fibers to motor units; 
neighbors_avoided = 3;
mf_assignment = assign_mf2mn(mf_centers, mu_centers, innerv_area, neighbors_avoided, Rmuscle);



%% Voronoi plot
%figure;
%plot_colored_voronoi2d(mf_centers, mf_assignment, Rmuscle);



%% Plot numbers of fibers per Motor Unit
figure;
histogram(mf_assignment,N); hold on; 
plot(innerv_num, 'r', 'linewidth', 1.5); 
title('Number of fibers distribution accross the units and target distribution');
xlim([0.5,N+0.5]); xlabel('Motor neuron'); ylabel('Number of fibers');
legend('Generated', 'Target');



%% Calculate and plot true innervation areas of motor neurons
mu_centers_result = zeros(N,2);
mu_areas_result = zeros(N,1);

phi_circle = linspace(0, 2*pi, 1000)';
phi_circle = phi_circle(1:end-1);
muscle_border = [Rmuscle * cos(phi_circle), Rmuscle * sin(phi_circle)];
figure; 
plot(muscle_border(:,1), muscle_border(:,2), 'k', 'linewidth', 1.2); hold on;

for i = 1:N
    mu_centers_result(i,:) = mean(mf_centers(mf_assignment == i, :));
    mu_areas_result(i) = sqrt(det(cov(mf_centers(mf_assignment == i, :))));
    
    %text(mu_centers_result(i,1), mu_centers_result(i,2), num2str(i));
    plot(mu_centers_result(i,1), mu_centers_result(i,2), 'b*');
    
    rad = mu_areas_result(i) * sqrt(chi2inv(.90, 2));
    mu_area_circle = [rad * cos(phi_circle) + mu_centers_result(i,1), ...
                      rad * sin(phi_circle) + mu_centers_result(i,2)];    
    plot(mu_area_circle(:,1), mu_area_circle(:,2), 'b');
    
    plot(mf_centers(mf_assignment == i,1), mf_centers(mf_assignment == i, 2), 'k.', 'MarkerSize', 5);
end

axis equal
title('Motor neuron innervation centers and areas over the muscle cross-section');
xlabel('x, mm'); ylabel('y, mm');
legend('Muscle border', 'Innervation centers', 'Innervation areas', 'Muscle fibers');



%% 
mu_areas_polygons = zeros(N,1);

figure; 
plot(muscle_border(:,1), muscle_border(:,2), 'k', 'linewidth', 1.2); hold on;

for i = 33:33
    assigned = find(mf_assignment == i);
    hp = plot(mf_centers(assigned,1), mf_centers(assigned, 2), '.', 'MarkerSize', 5);
    
    hull = convhull(mf_centers(assigned,1), mf_centers(assigned, 2));
    mu_areas_polygons(i) = polyarea(mf_centers(assigned(hull),1), mf_centers(assigned(hull),2));
    plot(mf_centers(assigned(hull),1), mf_centers(assigned(hull),2), '-', 'Color', hp.Color, 'linewidth', 1.2);
end

axis equal
title('Motor neuron innervation centers and areas over the muscle cross-section');
xlabel('x, mm'); ylabel('y, mm');
legend('Muscle border', 'Innervation centers', 'Innervation areas', 'Muscle fibers');

figure; 
plot(innerv_area); 
ylabel('Target innervation areas of motor neurons'); xlabel('Motor neuron');

yyaxis right; 
%plot((mu_areas_polygons-mean(mu_areas_polygons))/std(mu_areas_polygons)*std(innerv_area) + mean(innerv_area)); 
plot(mu_areas_polygons); hold on;
plot(mu_areas_result/std(mu_areas_result)*std(mu_areas_polygons), 'g-');
ylabel('Resulting innervation areas of motor neurons');
legend('Target areas of motor units', 'Result - convex hull areas', 'Result - SD of gaussian fit');




clear phi_circle rad mu_area_circle 




%% Theoretical diameters of muscle fibers (distribution parameters)
[diam_means, diam_stds] = get_mf_dist_parameters(innerv_area);


%% Assign model diameters to fibers
mf_diameters = zeros(size(mf_assignment));
for i = 1:N
    fibers_in_mu = find(mf_assignment == i);
    mf_diameters(fibers_in_mu) = diam_means(i) + diam_stds(i)*randn(size(fibers_in_mu));
end
relative_diameters = mf_diameters/mean(mf_diameters);


%% Assigning conduction velocities to fibers
mf_cv = assign_mf_cv(mf_diameters);
figure; histogram(mf_cv); title('Global distribution of MF conduction velocities');



%% Practical diameter distribution vs. model one
figure; 
%mf_diameters = 1e6*mf_diameters(sqrt(mf_centers(:,1).^2 + mf_centers(:,2).^2) < Rmuscle-0.1);
histogram(mf_diameters);
xlabel('Fiber diameters, $\mu$m');
ylabel('Historgram of fiber diameters over the muscle');
[muhat, sigmahat] = normfit(mf_diameters);
hold on; 
normfitx = linspace(min(mf_diameters), max(mf_diameters), 100);
plot(normfitx, sqrt(2*pi) * sigmahat * max(histcounts(mf_diameters)) * normpdf(normfitx, muhat, sigmahat), 'linewidth', 1.5);
text(60, 600, ['$\mu$ = ', num2str(muhat)], 'fontsize', 20);
text(60, 550, ['$\sigma$ = ', num2str(sigmahat)], 'fontsize', 20);
%title(sprintf('Dist of mf diameters, mean = %2.2fum, std = %2.2fum ', mean(mf_diameters)*1e6, std(mf_diameters)*1e6));

clear normfitx muhat sigmahat



%% Mean and std diameters vs motor unit indexes
diam_means_result = zeros(N,1);
diam_stds_result = zeros(N,1);
for i = 1:N
    diam_means_result(i) = mean(mf_diameters(mf_assignment == i));
    diam_stds_result(i) = std(mf_diameters(mf_assignment == i));
end
figure;
plot(diam_means_result * 1e6, 'linewidth',1.5); hold on;
plot(diam_means * 1e6, 'r', 'linewidth', 1.5)
xlabel('Motor unit'); ylabel('Mean fiber diameter, $\mu$m');
legend('Generated', 'Model');

figure; 
plot(diam_stds_result, 'linewidth', 1.5); hold on;
plot(diam_stds, 'linewidth', 1.5);
xlabel('Motor unit'); ylabel('Standard deviation of fiber diameters, $\mu$m');
legend('Generated', 'Model');

clear i



