
N = 100; %Number of mus
rr = 50; %Recruitment range: largest/smallest
rm = 0.75; %Recruitment maximum (when all the MUs are active)
Rmuscle = 5; %[mm]

mn_pool = MN_Pool_Sim(N, rr, rm);
mn_pool.distribute_innervation_centers(Rmuscle);

figure; ax = axes;
mn_pool.show_sizes(ax);

figure; ax = axes;
mn_pool.show_centers(ax, Rmuscle);

%% Define EXC-FR curves
mn_pool.set_deluca_mdl('fdi');
%mn_pool.generate_minfr('linear_rt', -10, 14); % Parameter 1: slope, Parameter 2: intersect
%mn_pool.generate_maxfr('linear_rt', -10, 25);
%mn_pool.generate_frs('linear_rt', -25, 50);

mn_pool.show_fr_exc_curves();

mn_pool.CV = 1/6;

clear rr rm