function [ diameters ] = estimate_real_diameters( mf_centers, nn_or_border, mode)
% function [ diameters ] = estimate_real_diameters( mf_centers, par1, mode)
% This function receives the muscle fibers' centers and estimates an actual
% diameter of each fiber as function of it's position relative to its
% neighbours. 
% Two modes: 
% 'simple' provides faster result, diameter is an average distance to par1
% nearest neighbours. 
% 'voronoi' calculates equivalent diameter as:
% d = sqrt(4/pi * voronoi_polygone_area); 
% a bit slower, but more reliable, though produces some outliers due
% to border effect.



if nargin < 3
    mode = 'voronoi';
end

% Going from millimeters to meters
mf_centers = mf_centers/1000;

switch mode
    case 'simple'
        n_neighbours = nn_or_border;
        disp('Simpler actual diameter estimation through dist to nearest neighbours');
        distmat = pdist2(centers,centers);
        distmat(distmat == 0) = Inf;
        distmat_sorted = sort(distmat, 'ascend');
        diameters = mean(distmat_sorted(1:n_neighbours, :));
    
    case 'voronoi'
        muscle_border = nn_or_border;
        muscle_border = muscle_border/1000;
        nmf = length(mf_centers);
        diameters = zeros(nmf,1);
        
        % Get Voronoi diagram of the muscle fibers distribution
        % Add muscle border points to prevent infinite areas for the Voronoi cells
        % closest to the border:
        disp('Generating the voronoi map of muscle fiber distribution...');
        [v,c] = voronoin([mf_centers; muscle_border], {'Qbb'});
        c = c(1:nmf);
        
        for i = 1:nmf
            polygon_vertices = v(c{i},:);
            voronoi_cell_area = polyarea(polygon_vertices(:,1), polygon_vertices(:,2));
            diameters(i) = sqrt(4/pi * voronoi_cell_area);
        end
end



