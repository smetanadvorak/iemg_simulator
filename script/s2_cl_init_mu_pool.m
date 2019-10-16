%% Glossary:
% generate - generate data/parameters/entities according to a method
% defined elsewhere.
% init - kind of generate, but if its called the first for a new object.
% calc - calculate new data as a transformation of a previously generated data.

%% 
mu_pool = MU_Pool_Sim(mn_pool); % Initialize the MU_Pool object with motor neuron pool
clear mn_pool

%% Generate MFs
mu_pool.generate_mfs(Rmuscle, Dmf); 
mu_pool.show_mf_centers();

%% Calculate theoretical innervation areas for MNs
overlap_degree = mu_pool.calc_innervation_areas();

%% Calculate theoretical innervation numbers for MNs (number of fibers innervated by each MN)
mu_pool.calc_innervation_numbers();

%% Assign MFs to MNs
mu_pool.assign_mfs2mns(4); % Greater the parameter is, less the fibers of the same MU are allowed to be adjacent to each other. See the (Akhmadeev et al. 2019) paper, adjacency parameter.

%% Show results of innervation
mu_pool.calc_innervation_numbers_res(); 

mu_pool.calc_innervation_areas_res('polygone_area', 0.95); % Parameter defines the way the resulting innervation area is calculated.
                                                           %'root_variance'; %'confidence_ellipse'; 'polygone_area'; 

mu_pool.show_innervation_numbers();
mu_pool.show_innervation_areas_1d();
mu_pool.show_innervation_areas_2d(1:2:N); 
title('');

%% Generate diameters and conduction velocities
mu_pool.generate_mf_diameters();
mu_pool.generate_mf_cvs();
mu_pool.show_diameters_distribution();

%%

clear Dmf Nmf
