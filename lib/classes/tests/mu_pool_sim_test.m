%% Muscle geometry
addpath(genpath('../..'));

% Everything is in millimeters. Muscle in vicinity of the electrode is approximated as a cylinder with
Lmuscle = 40; % [mm]
Dmf = 400; % Density of muscle fibers per square millimetre (Hamilton-Wright 2005)

N = 120; %Number of mus
rr = 50; %Recruitment range: largest/smallest
rm = 0.75; %Recruitment maximum (when all the MUs are active)
Rmuscle = 5; %[mm]

%% First, create an MN pool object
mn_pool = MN_Pool_Sim(N, rr, rm);
mn_pool.distribute_innervation_centers(Rmuscle);

%% 
mu_pool = MU_Pool_Sim(mn_pool);

%% Generate MFs
mu_pool.generate_mfs(Rmuscle, Dmf);
mu_pool.show_mf_centers();

%% Innervation Area
overlap_degree = mu_pool.calc_innervation_areas();

%% Innervation Numbers
mu_pool.calc_innervation_numbers();
mu_pool.calc_innervation_areas();

%% Assign mfs to mns
mu_pool.assign_mfs2mns(3);

%% Show results
mu_pool.calc_innervation_numbers_res();

ia_type = 'polygone_area';%'root_variance'; %'confidence_ellipse'; ; 
mu_pool.calc_innervation_areas_res(ia_type);

mu_pool.show_innervation_numbers();
mu_pool.show_innervation_areas_1d();
mu_pool.show_innervation_areas_2d();
