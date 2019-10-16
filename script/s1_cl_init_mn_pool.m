
%% Muscle geometry
% Everything is in millimeters. Muscle in vicinity of the electrode is approximated as a cylinder with
Dmf = 400; % Density of muscle fibers per square millimetre (Hamilton-Wright 2005)
Nmf = 34000; % Expected number of muscle fibers in the muscle (40k for FDI, see Feinstein - Morphologic studies ... 1995)
Lmuscle = 30; % [mm]
Rmuscle = 5;%sqrt((Nmf/Dmf)/pi); %[mm]

%% MN pool parameters
N = 100; %Number of mus (120 for FDI, see Feinstein - Morphologic studies ... 1995)
rr = 50; %Magnitude of RT distribution: largest/smallest
rm = 0.75; %Recruitment maximum (when all the MUs are active)

mn_pool = MN_Pool_Sim(N, rr, rm);
mn_pool.distribute_innervation_centers(Rmuscle);

figure; ax = axes;
mn_pool.show_sizes(ax);

figure; ax = axes;
mn_pool.show_centers(ax, Rmuscle);

%% Define EXC-FR curves
%mn_pool.set_deluca_mdl('fdi');
mn_pool.generate_minfr('linear_rt', -5, 10); % Parameter 1: slope, Parameter 2: intersect
mn_pool.generate_maxfr('linear_rt', -10, 40);
mn_pool.generate_frs('linear_rt', -20, 50);

mn_pool.show_fr_exc_curves();

mn_pool.CV = 1/7;

clear rr rm