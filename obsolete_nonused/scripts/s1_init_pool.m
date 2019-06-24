

N = 100; % Number of mus
rr = 50; % Recruitment range: largest/smallest
rm = 0.75; % Recruitment maximum (when all the MUs are active)
Rmuscle = 5; % [mm]s


%% MU Recruitment threshold
[sz, rt] = generate_mu_recruitment_thresholds(N, rr, rm);
figure; stem(1:N, sz); hold on; 

title('Thresholds and sizes distribution over the motor neuron pool');
xlabel('Index of MU');
ylabel('Motor neuron size, normalized units');
yyaxis right;
stem(1:N,rt, 'o');
ylabel('Excitation rate, normalized units');
legend({'MU sizes','Recruitment thresholds'});
yyaxis left;
figure2page;
% Check: smallest/largest in terms of mf number should be ~ 1/100
% (Fuglevand, p.2476)


%% MU Firing rate models
[minfr] = generate_mu_minfr(sz, 'constant', 8);
[maxfr] = generate_mu_maxfr(sz, 'linear_rt', 20, 20);
%[frs] = generate_mu_fr_slope(sz, 'constant');
frs = (maxfr - minfr)./(1-rt); % Maximal excitation is 1, maxfr is occurring only at e=1 (100% of MVC)


%% MU innervation area
innerv_area = sz/sum(sz) * pi * Rmuscle^2;


%% Generate innervation centres
% Uniform distribution of MU centres.
% Coordinates are polar and centered at muscle cross-section center.
% Generated in cartesian
% no_mu_land = 0; %sqrt(mean(innerv_area)/pi)/4;
% fp_rad = 0.5;
% cd ./voronoi_sampling/bogdan_voronoi/build
% system('rm voronoi'); system('make');
% command = sprintf('./voronoi %d %1.2f', N, fp_rad);
% system(command);
% mu_centers = csvread('voronoi.csv');
% system('rm voronoi.csv');
% cd ../../..
% mu_centers = (mu_centers - 5);
% mu_centers = mu_centers(sqrt(mu_centers(:,1).^2 + mu_centers(:,2).^2) < fp_rad, :);
% mu_centers = mu_centers/fp_rad * (Rmuscle - no_mu_land);
% mu_centers = mu_centers(end:-1:1, :);
% 
% clear no_mu_land fp_rad command
%cd geodesic_farthest_point
%addpath('toolbox_general/');
%addpath('toolbox_graph/');
%addpath('toolbox_signal/');

n = 256;
density_map = ones(n);
[X,Y] = meshgrid(1:n,1:n);
density_map(sqrt((X-n/2).^2 + (Y-n/2).^2) > n/2-1) = 0;

vertices = zeros(2,N);
vertices(:,1) = [1;1];
for i = 2:(N+1)
    D = perform_fast_marching(1./density_map, vertices(:, 1:i-1));
    [~,ind] = max(D(:));
    [x,y] = ind2sub([n n],ind);
    vertices(:,i) = [x;y];
end

mu_centers = vertices(:,end:-1:2)';
mu_centers = (mu_centers - n/2) / max(sqrt((mu_centers(:,1)-n/2).^2 + (mu_centers(:,2)-n/2).^2)) * Rmuscle;


%% Plot innervation centers and areas
clear n density_map X Y vertices

phi_circle = linspace(0, 2*pi, 1000)';
phi_circle = phi_circle(1:end-1); % No overlapping
muscle_border = [Rmuscle * cos(phi_circle), Rmuscle * sin(phi_circle)];
figure; 
plot(muscle_border(:,1), muscle_border(:,2), 'b', 'linewidth', 1.2); hold on;
for i = 1:N
    text(mu_centers(i,1), mu_centers(i,2), num2str(i));
    rad = sqrt(innerv_area(i)/pi);
    mu_area_circle = [rad * cos(phi_circle) + mu_centers(i,1), ...
                      rad * sin(phi_circle) + mu_centers(i,2)];
    
    plot(mu_area_circle(:,1), mu_area_circle(:,2), 'k');
end
axis equal
%title('Motor neuron innervation centers and areas over the muscle cross-section');
xlabel('x, mm'); ylabel('y, mm');
legend('Muscle border', 'Innervation areas');

clear phi_circle rad mu_area_circle


