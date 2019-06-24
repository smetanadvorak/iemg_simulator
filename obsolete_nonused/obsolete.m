%% Old stuff:
%% Muscle fibers coordinates generation
% 
% % =========== Fiber diameters due to MU's size modification ===========
% % Theoretical diameters of muscle fibers (distribution parameters)
% [diam_means, diam_stds] = get_mf_dist_parameters(innerv_area);
% 
% voronoi_sampling_border_rad = 4; % Voronoi sampling is originally done over a square domain (10x10). This is the radius of a circular zone in the middle of the domain. Smaller radii negate the border effects but also require longer computation times.
% 
% cd ./voronoi_sampling/bogdan_voronoi
% %csvwrite('mu_centers.csv', mu_centers);
% %csvwrite('mu_areas.csv', innerv_area);
% %csvwrite('mu_normalized_mean_diameters.csv', diam_means/55e-6);
% 
% cd ./build
% system('rm voronoi'); system('make');
% 
% command = sprintf('./voronoi %d %1.2f', Nmf, voronoi_sampling_border_rad); system(command);
% mf_centers = csvread('./voronoi.csv'); system('rm voronoi.csv');
% system('rm *.csv'); cd ..
% system('rm *.csv');
% cd ../..
% 
% mf_centers = (mf_centers - 5);
% mf_centers = mf_centers(sqrt(mf_centers(:,1).^2 + mf_centers(:,2).^2) < voronoi_sampling_border_rad, :);
% mf_centers = mf_centers/voronoi_sampling_border_rad * Rmuscle;
% clear voronoi_sampling_border_rad command 



%% Plot MF diameters distribution
% figure; ax =[]; set(gcf, 'Name', 'MF Diameters Distribution, um');
% side = ceil(sqrt(N));
% for i = 1:side
%     for j = 1:side
%         ind = (i-1)*side + j;
%         ax(ind) = subplot(side, side, ind);
%         if ind <= N
%             histogram(mf_diameters(mf_assignment == ind) * 1e6)
%             axis tight; 
%             title(sprintf('MD: %2.2f', mean(mf_diameters(mf_assignment == ind)) * 1e6));
%         end
%     end
% end
% clear side i j ind ax


%% Calculating actual diameters of fibers (under question)
% phi_circle = linspace(0, 2*pi, 10000)'; %Should be sampled as finely as possible;
% phi_circle = phi_circle(1:end-1); % No overlapping
% muscle_border = [Rmuscle * cos(phi_circle), Rmuscle * sin(phi_circle)];
% mf_diameters = estimate_real_diameters(mf_centers, muscle_border, 'voronoi');
% 
% % Redistribute the mf_diameters of fibers on the border of the muscle,
% % their diameters are way underestimated:
% diam_mean_global = mean(mf_diameters);
% diam_std_global = std(mf_diameters);
% mf_diameters(mf_diameters < diam_mean_global - 3*diam_std_global) = diam_mean_global + diam_std_global*randn;
% %relative_actual_mf_diameters = mf_diameters/mean(mf_diameters);
% clear phi_circle