%% Glossary:
% generate - generate data/parameters/entities according to a method
% defined elsewhere.
% init - kind of generate, but if its called the first for a new object.
% calc - calculate new data as a transformation of a previously generated data.



%% Muscle geometry
% Everything is in millimeters. Muscle in vicinity of the electrode is approximated as a cylinder with
Lmuscle = 50; % [mm]
Dmf = 400; % Density of muscle fibers per square millimetre (Hamilton-Wright 2005)
Nmf = round((Rmuscle^2) * pi * Dmf); % Expected number of muscle fibers in the muscle. 

%% 
mu_pool = MU_Pool_Sim(mn_pool);
clear mn_pool

%% Generate MFs
mu_pool.generate_mfs(Rmuscle, Dmf);
mu_pool.show_mf_centers();

%% Innervation Area
overlap_degree = mu_pool.calc_innervation_areas();

%% Innervation Numbers
mu_pool.calc_innervation_numbers();

%% Assign mfs to mns
mu_pool.assign_mfs2mns(5);


%% Show results for innervation
mu_pool.calc_innervation_numbers_res();

mu_pool.calc_innervation_areas_res('polygone_area', 0.95); %'root_variance'; %'confidence_ellipse'; 'polygone_area'; 

mu_pool.show_innervation_numbers();
mu_pool.show_innervation_areas_1d();
mu_pool.show_innervation_areas_2d(randperm(mu_pool.N, 50));
title('');

%% Generate diameters and conduction velocities
mu_pool.generate_mf_diameters();
mu_pool.generate_mf_cvs();
mu_pool.show_diameters_distribution();

%%

clear Dmf Nmf
