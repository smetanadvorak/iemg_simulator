
open innervation_numbers.fig
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', [0, 0, 20, 15]); %[0, 0, 7.4, 5.2], [0, 0, 10.5, 7.4]);
set(gcf, 'papersize', [21 29.7]);
set(gca, 'fontunits', 'points', 'fontsize', 20);

%print -dpng -r600 -painters mn_subsets.png
print -depsc2 -painters n_fibers.eps

%% 
open territories.fig
axis tight
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', [0, 0, 20, 15]); %[0, 0, 7.4, 5.2], [0, 0, 10.5, 7.4]);
set(gcf, 'papersize', [21 29.7]);
set(gca, 'fontunits', 'points', 'fontsize', 20);

%print -dpng -r600 -painters mn_subsets.png
print -depsc2 -painters territories.eps

%% 
open mn_centers.fig
axis tight
set(gca, 'xtick', [-4:2:4]); set(gca, 'ytick', [-5:1:5]);

set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', [0, 0, 20, 15]); %[0, 0, 7.4, 5.2], [0, 0, 10.5, 7.4]);
set(gcf, 'papersize', [21 29.7]);
set(gca, 'fontunits', 'points', 'fontsize', 18);

%print -dpng -r600 -painters mn_subsets.png
print -depsc2 -painters mn_centers.eps