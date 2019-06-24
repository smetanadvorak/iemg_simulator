N = 100; %Number of mus
rr = 50; %Recruitment range: largest/smallest
rm = 0.70; %Recruitment maximum (when all the MUs are active) %0.67 for FDI from deLuca
Rmuscle = 5; %[mm]

mn_pool = MN_Pool_Sim(N, rr, rm);
mn_pool.distribute_innervation_centers(Rmuscle);

figure; ax = axes;
mn_pool.show_sizes(ax);

figure; ax = axes;
mn_pool.show_centers(ax, Rmuscle);

%%
mn_pool.set_deluca_mdl('fdi');

mn_pool.show_fr_exc_curves();
