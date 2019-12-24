classdef MU_Pool_Sim < handle
    %MU_POOL_SIM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N
        mn_pool;
        innervation_areas;
        innervation_areas_res;
        
        muscle_area2max_ia;
        
        innervation_numbers;
        innervation_numbers_res;
        
        mf_centers;
        Nmf;
        Dmf; 
        Rmuscle;
        muscle_border;
        
        mf_diameters
        mf_cv
        
        assignment;
    end
    
    methods
        function obj = MU_Pool_Sim(mn_pool)
            obj.mn_pool = mn_pool;
            obj.N = mn_pool.N;
        end
        
        function generate_mfs(obj, Rmuscle, Dmf)
            if nargin < 3
                Dmf = 400; % Density of muscle fibers per square millimetre (Hamilton-Wright 2005)
            end
            
            obj.Dmf = Dmf; obj.Rmuscle = Rmuscle;
            
            % Expected number of muscle fibers in the muscle. 
            obj.Nmf = round((Rmuscle^2) * pi * Dmf); 

            
            % Get the class file path, the precalculated data is there;
            cls_path = fileparts(mfilename('fullpath'));
            cur_path = pwd;
            
            % Download the precalculated distribution
            cd(cls_path)
                obj.mf_centers = csvread('./mf_data/voronoi_pi1e5.csv');
            cd(cur_path)
            
            % Adjust the loaded centers to the expected number of fibers
            % and to muscle radius
            obj.mf_centers = (obj.mf_centers - 5)/4; % 4 may be unnecessary here
            [dists, inds] = sort(sqrt(obj.mf_centers(:,1).^2 + obj.mf_centers(:,2).^2));
            obj.mf_centers = obj.mf_centers(inds(1:obj.Nmf), :) / dists(obj.Nmf+1) * Rmuscle;
            
            % Create muscle border for plotting
            phi_circle = linspace(0, 2*pi, 1000)';
            phi_circle = phi_circle(1:end-1); 
            obj.muscle_border = [obj.Rmuscle * cos(phi_circle), obj.Rmuscle * sin(phi_circle)];
        end
        
        function overlap_degree = calc_innervation_areas(obj, muscle_area2max_ia)
            if(nargin < 2)
                muscle_area2max_ia = 4;
            end
            obj.muscle_area2max_ia = muscle_area2max_ia;
            
            muscle_area = pi * obj.Rmuscle^2;
            % Innervation areas are proportional to the sizes of MNs 
            % Areas are scaled in such way that the larges MN innervate 
            % 1/max_ia2muscle_area of the total muscle area.
            obj.innervation_areas = obj.mn_pool.sz/max(obj.mn_pool.sz) * muscle_area / muscle_area2max_ia;
            overlap_degree = sum(obj.innervation_areas)/muscle_area;
        end
        
        function calc_innervation_numbers(obj)
            muscle_area = pi * obj.Rmuscle^2;
            obj.innervation_numbers = round(obj.innervation_areas/sum(obj.innervation_areas) * muscle_area * obj.Dmf);
        end
        
        function assign_mfs2mns(obj, n_neighbours, conf)
            % This function assigns fibers to neurons according to several rules:
            % 1) Proximity to the innervation center;
            % 2) Innervation area and target number of innervated fibers)
            % 3) Self-avoiding phenomena (see n_neighbours variable)
            
            if nargin < 2
                n_neighbours = 3;
            end
            
            % Confidence interval defines the ratio between the innervation
            % area and the variance of the gaussian distribution. The
            % larger, the tighter the clusters are.
            if nargin < 3
                conf = 0.999;
            end
            
            % Out-of-muscle area compensation
            % Calculates how much of the MU's gaussian distribution is outside of the
            % muscle border and inflates the rest of the distribution according to it
            borderfun_pos = @(x)real( sqrt(obj.Rmuscle^2 - x.^2));
            borderfun_neg = @(x)real(-sqrt(obj.Rmuscle^2 - x.^2));
            out_circle_coeff = ones(obj.N,1);
            
            c =  chi2inv(conf,2);
            sigma = @(ia) eye(2) * ia / pi / c;
            
            for mu = 1:obj.N
                probfun = @(x,y) reshape(mvnpdf([x(:), y(:)], obj.mn_pool.centers(mu,:), sigma(obj.innervation_areas(mu))), size(x)); % adapted to integral2
                in_circle_int = integral2(probfun, -obj.Rmuscle, obj.Rmuscle, borderfun_neg, borderfun_pos);
                out_circle_coeff(mu) = 1/in_circle_int;
            end

            % Distmat for supressing closest neighbours
            if n_neighbours
                neighbours = knnsearch(obj.mf_centers, obj.mf_centers, 'K', n_neighbours+1);
                neighbours = neighbours(:,2:end);
            end
            
            % Assignment procedure
            obj.assignment = nan(obj.Nmf,1);
            randomized_mf = randperm(obj.Nmf);
            
            i = 0;
            for mf = randomized_mf
                probs = zeros(obj.N,1);
                for mu = 1:obj.N
                    % Supression assignment if neighbours are from the same MU
                    % Promotes intermingling
                    % Problem: spreads the MUs too much;
                    if n_neighbours && any(obj.assignment(neighbours(mf,:)) == mu)
                        probs(mu) = 0;
                    else
                        % A priori probability of the assignment
                        apriori_prob = obj.innervation_numbers(mu)/obj.Nmf;
                        
                        % Likelihood coming from clustered nature of mf distribution
                        clust_hood = mvnpdf(obj.mf_centers(mf,:)', obj.mn_pool.centers(mu,:)', sigma(obj.innervation_areas(mu)));
                        clust_hood = clust_hood * out_circle_coeff(mu);
                        
                        % Final a posteriori probability
                        probs(mu) = apriori_prob * clust_hood;
                    end
                end
                
                probs = probs/sum(probs);
                obj.assignment(mf) = randsample(obj.N, 1, true, probs);
                
                 i = i+1;
                if ~mod(i,1000)
                    fprintf('%d of muscle fibers assigned\n', i);
                end
            end
            
            obj.calc_innervation_numbers_res();
            obj.calc_innervation_areas_res();
        end
        
        function innervation_numbers_res = calc_innervation_numbers_res(obj)
            obj.innervation_numbers_res = histcounts(obj.assignment, obj.N);
            innervation_numbers_res = obj.innervation_numbers_res;
        end
        
        function innervation_areas_res = calc_innervation_areas_res(obj, type, conf)
            if nargin < 2
                type = 'confidence_ellipse';
            end
            if nargin < 3 && strcmp(type, 'confidence_ellipse')
                conf = 0.95;
            end
            
            obj.innervation_areas_res = zeros(obj.N,1);
            
            switch type
                case 'confidence_ellipse'
                    for m = 1:obj.N
                        covariance = cov(obj.mf_centers(obj.assignment == m, :));
                        [~, eigenval] = eig(covariance);
                        
                        chisquare_val = chi2inv(conf, 2);
                        % Innervation area is area of the confidence interval
                        % ellipse
                        obj.innervation_areas_res(m) = prod(sqrt(diag(eigenval))) * chisquare_val * pi;
                    end
                    
%                 case 'root_variance'
%                     for m=1:obj.N
%                         %[~, sigma] = normfit(obj.mf_centers(obj.assignment == m, :)));
%                         %obj.innervation_areas_res(m) = sqrt(cummult(diag(sigma)));
%                         covariance = cov(obj.mf_centers(obj.assignment == m, :));
%                         obj.innervation_areas_res(m) = sqrt(prod(diag(eig(covariance))));
%                     end
                    
                case 'polygone_area'
                    for m = 1:obj.N
                        assgn = find(obj.assignment == m);
                        hull = convhull(obj.mf_centers(assgn,1), obj.mf_centers(assgn, 2));
                        obj.innervation_areas_res(m) = polyarea(obj.mf_centers(assgn(hull),1), obj.mf_centers(assgn(hull),2));
                    end
            end
            
            innervation_areas_res = obj.innervation_areas_res;
        end
        
        %% Diameters and coducton velocities methods
        
        function generate_mf_diameters(obj)
            [diam_means, diam_stds] = get_mf_dist_parameters(obj.mn_pool.sz);
            %% Theoretical diameters of muscle fibers (distribution parameters)
            
            %% Assign model diameters to fibers
            % A method by Stashuk and Hamilton-Wright (PHYSIOLOGICALLY BASED SIMULATION OF CLINICAL EMG SIGNALS)
            obj.mf_diameters = zeros(size(obj.assignment));
            for m = 1:obj.N
                fibers_in_mu = find(obj.assignment == m);
                obj.mf_diameters(fibers_in_mu) = diam_means(m) + diam_stds(m)*randn(size(fibers_in_mu));
            end
        end
        
        function generate_mf_cvs(obj)
            %% Assign conduction velocities to fibers
            obj.mf_cv = assign_mf_cv(obj.mf_diameters);
        end
        
        
        %% Visualization methods
        function show_mf_centers(obj, ax)
            if nargin < 2
                figure;
                ax = axes();
            end
            % Draw circle for muscle border
            phi_circle = linspace(0, 2*pi, 1000)';
            phi_circle = phi_circle(1:end-1); 
            muscle_border = [obj.Rmuscle * cos(phi_circle), obj.Rmuscle * sin(phi_circle)];
            
            plot(ax, muscle_border(:,1), muscle_border(:,2), 'k'); hold on;
            plot(ax, obj.mf_centers(:,1), obj.mf_centers(:,2), 'k.', 'markersize',1); axis equal;
            xlabel('X, mm'); ylabel('Y, mm'); title('Muscle fibers'' centers in muscle cross-section');
        end
        
        function show_innervation_areas_1d(obj, ax)
            if nargin < 2
                figure;
                ax = axes();
            end
            plot(obj.innervation_areas, 'r', 'linewidth', 1.5); hold on;
            bar(obj.innervation_areas_res, 'b');
            %title('Innervation areas of the units and their target values');
            xlim([0.5,obj.N+0.5]); xlabel('Motor neuron'); ylabel('Innervation areas, $mm^2$');
            legend('Target', 'Generated', 'location', 'nw');            
        end
        
        function hps = show_innervation_areas_2d(obj, inds, ax)
            % Show innervation areas of chosen motor units
            if nargin < 2
                inds = round(linspace(1,obj.N,10));
            end
            if nargin < 3
                figure; ax = axes;
            end
            axis(ax);
            
            % Draw circle for muscle border
            plot(ax, obj.muscle_border(:,1), obj.muscle_border(:,2), 'k', 'linewidth', 1); hold on;
            hps = [];
            for m = 1:numel(inds)
                assgn = find(obj.assignment == inds(m));
                hull = convhull(obj.mf_centers(assgn,1), obj.mf_centers(assgn, 2));
                hp = plot(obj.mf_centers(assgn(hull),1), obj.mf_centers(assgn(hull),2), '-', 'linewidth', 0.75); hold on;
                plot(obj.mf_centers(assgn,1), obj.mf_centers(assgn, 2), '.', 'MarkerSize', 0.1, 'Color', hp.Color); 
                hps = [hps, hp]; 
            end
            for m = 1:numel(inds)
                %text(obj.mn_pool.centers(inds(m),1), obj.mn_pool.centers(inds(m),2), num2str(inds(m)), 'Color', hps(m).Color, 'fontsize', 7, 'linewidth',0.5, 'horizontalalignment', 'center');
            end
            
            axis equal
            %title('Motor neuron innervation centers and areas over the muscle cross-section');
            xlabel('x, mm'); ylabel('y, mm');
            %legend('Muscle border', 'Innervation centers', 'Innervation areas', 'Muscle fibers');
            
        end
       
        
        function show_innervation_numbers(obj, ax)
            if nargin < 2
                figure;
                ax = axes();
            end
            
            histogram(ax, obj.assignment, obj.N); hold on;
            plot(obj.innervation_numbers, 'r', 'linewidth', 1.5);
            %title('Number of fibers distribution accross the units and the target distribution');
            xlim([0.5,obj.N+0.5]); xlabel('Motor neuron'); ylabel('Number of fibers');
            legend('Generated', 'Target', 'location', 'nw');
        end
        
        function show_cv_distribution(obj, ax)
            if nargin < 2
                figure;
                ax = axes();
            end
            histogram(ax, obj.mf_cv); 
            %title('Global distribution of MF conduction velocities');
        end
        
        function show_diameters_distribution(obj, ax)
            if nargin < 2
                figure;
                ax = axes();
            end
            histogram(ax, obj.mf_diameters * 10^6); 
            xlabel('Fiber diameters, $\mu$m');
            ylabel('Historgram of fiber diameters over the muscle');
            [muhat, sigmahat] = normfit(obj.mf_diameters * 10^6);
            hold on;
            normfitx = linspace(min(obj.mf_diameters * 10^6), max(obj.mf_diameters * 10^6), 100);
            %plot(normfitx, sqrt(2*pi) * sigmahat * max(histcounts(obj.mf_diameters * 10^6)) * normpdf(normfitx, muhat, sigmahat), 'linewidth', 1.5);
            plot(normfitx, sqrt(2*pi) * 9 * max(histcounts(obj.mf_diameters * 10^6)) * normpdf(normfitx, 55, 9), 'linewidth', 2);
            text(70, 900, ['$\hat{\mu}$ = ', num2str(muhat)], 'fontsize', 20);
            text(70, 800, ['$\hat{\sigma}$ = ', num2str(sigmahat)], 'fontsize', 20);
            legend('Resulting distribution', 'Experimental distribution');
        end
    end
end

