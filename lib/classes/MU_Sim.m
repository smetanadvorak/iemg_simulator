classdef MU_Sim < handle
    
    properties
        
        muap;
        sfaps;
        
        % Muscle fibers properties
        Nmf;
        mf_centers; % muscle fibers' x,y positions;
        mf_left_end; % coordinates of muscle fibers left ends;
        mf_right_end; % coordinates of muscle fibers right ends;
        mf_diameters; % diameters of muscle fibers;
        mf_cv; % conduction velocities
        
        % Neuromuscular junctions
        nmj_z; % zs of neuromuscular junctions;
        branch_points_xy;
        branch_points_z;
        nmj_cv; % nerves' branches' conduction velocities (as many values as layers in the model);
        mnap_delays;
        nerve_paths;
        
        % Jitter
        nmj_jitter; %value of neuromuscular jitter
        
        % Innervation geometry
        actual_center;
        nominal_center;
        
        % Simulation parameters
        dt;
        dz;
        
        % Electrode parameters
        Npt; %Number of points
        
    end
    
    methods
        %% Constructor
        function obj = MU_Sim(mf_centers, mf_left_end, mf_right_end, mf_diameters, mf_cv, nmj_cv)
            if nargin ~= 0
                obj.Nmf = size(mf_centers,1);
                obj.mf_centers = mf_centers;
                obj.mf_left_end = obj.distributeParameter(mf_left_end);
                obj.mf_right_end = obj.distributeParameter(mf_right_end);
                obj.mf_diameters = obj.distributeParameter(mf_diameters);
                obj.mf_cv = obj.distributeParameter(mf_cv);
                obj.nmj_cv = nmj_cv;
                % Automatically initialized:
                %obj.nmj_z = (obj.mf_right_end - obj.mf_left_end)/2;
            end
        end
        
        %% SFAPs calculation function
        function calc_sfaps(obj, dt, dz, pts_world, pts_normals, min_radial_dist)
            if nargin < 6
                min_radial_dist = mean(obj.mf_diameters) * 1000; % Get millimeters
            end
            obj.Npt = size(pts_world, 1);
            
            obj.dt = dt;
            obj.dz = dz;
            
            % Get the grid that covers the whole length of fiber and most
            % of the action potential timespan;
            t = 0:dt:(2 * max([(obj.nmj_z - obj.mf_left_end)./obj.mf_cv ; (obj.mf_right_end - obj.nmj_z)./obj.mf_cv]));
            t = t';
            obj.sfaps = zeros(numel(t), obj.Npt, obj.Nmf);
            
            % Simulation loop
            for fb = 1:obj.Nmf
                % Get current density of selected fiber:
                % Symmetrical z relative to the current nmj coordinate (gives symmetrical spatial current densities):
                z_left = obj.nmj_z(fb):-dz:(obj.mf_left_end(fb));
                z_right = obj.nmj_z(fb):dz:(obj.mf_right_end(fb));
                z = transpose([z_left(end:-1:1), z_right(2:end)]);
                mf_coord_3d = [repmat(obj.mf_centers(fb,:), numel(z), 1), z(:)];
                
                % Get current density of selected fiber
                current_density = get_current_density_symm(t, z, obj.nmj_z(fb), obj.mf_right_end(fb) - obj.nmj_z(fb), obj.nmj_z(fb) - obj.mf_left_end(fb), obj.mf_cv(fb), obj.mf_diameters(fb));
                
                for pt = 1:obj.Npt
                    radial_dist = sqrt(sum((pts_world(pt,1:2) - obj.mf_centers(fb,:)).^2));
                    radial_dist(radial_dist < min_radial_dist) = min_radial_dist;
                    response_to_elem_current = get_elementary_current_response(z, pts_world(pt,3), radial_dist);
                    %if norm(pts_normals(pt,:)) > 0
                    %    effective_elem_area = get_visible_area(mf_coord_3d, pts_world(pt,:), pts_normals(pt,:));
                    %else
                        effective_elem_area = 1;
                    %end
                    obj.sfaps(:,pt,fb) = current_density' * (response_to_elem_current .* effective_elem_area);
                end
            end
            obj.shift_sfaps(dt);
        end
 
        
        %% Precalculated SFAPs
        function calc_sfaps_precalc(obj, dt, dz, pts_world, pts_normals)
            obj.Npt = size(pts_world, 1);
            
            obj.dt = dt;
            obj.dz = dz;
            
            % Get the grid that covers the whole length of fiber and most
            % of the action potential timespan;
            t = 0:dt:(2 * max([(obj.nmj_z - obj.mf_left_end)./obj.mf_cv ; (obj.mf_right_end - obj.nmj_z)./obj.mf_cv]));
            t = t';
            
            % Get current density of selected fiber:
            % Symmetrical z relative to the current nmj coordinate (gives symmetrical spatial current densities):
            z = transpose(0:dz:max(max(obj.mf_right_end), max(obj.mf_cv)*max(t)));

            precalculated_cd = get_tm_current_dz(0:dz:(10*max(z)));
            
            
            obj.sfaps = zeros(numel(t), obj.Npt, obj.Nmf);
            
            % Simulation loop
            for fb = 1:obj.Nmf                
                % Get current density of selected fiber
                mf_coord_3d = [repmat(obj.mf_centers(fb,:), numel(z), 1), z(:)];
                
                current_density = get_current_density_fast(precalculated_cd, t, z, obj.nmj_z(fb), obj.mf_right_end(fb) - obj.nmj_z(fb), obj.nmj_z(fb) - obj.mf_left_end(fb), obj.mf_cv(fb), obj.mf_diameters(fb));
                for pt = 1:obj.Npt
                    response_to_elem_current = get_electric_distance(z, pts_world(pt,3), sqrt(sum((pts_world(pt,1:2) - obj.mf_centers(fb,:)).^2)));
                    if norm(pts_normals(pt,:)) > 0
                        effective_elem_area = get_visible_area(mf_coord_3d, pts_world(pt,:), pts_normals(pt,:));
                    else
                        effective_elem_area = 1;
                    end
                    obj.sfaps(:,pt,fb) = current_density' * (response_to_elem_current .* effective_elem_area);
                end
            end
            obj.shift_sfaps(dt);
        end
 
        
        
        
        
        %%
        function shift_sfaps(obj, dt)
            obj.calc_mnap_delays;
            for fb = 1:obj.Nmf
                for pt = 1:obj.Npt
                    obj.sfaps(:,pt,fb) = shift_padding(obj.sfaps(:,pt,fb), floor(obj.mnap_delays(fb)/dt), 1);
                    obj.sfaps(:,pt,fb) = hr_shift_template(obj.sfaps(:,pt,fb), mod(obj.mnap_delays(fb), dt));
                end
            end
        end
        
        %% MUAP-related functions
        % ToDo: merge shift_sfaps with calc_muap
        function calc_mnap_delays(obj)
            obj.mnap_delays = obj.nerve_paths./repmat(obj.nmj_cv(:)', obj.Nmf, 1);
            obj.mnap_delays = sum(obj.mnap_delays,2);
        end
        
        function muap = calc_muap(obj, jitter_std)
            if nargin < 2
                jitter_std = obj.nmj_jitter;
            end
            
            if jitter_std ~= 0
                delays = jitter_std * randn(obj.Nmf,1);
                jittered_sfaps = zeros(size(obj.sfaps));
                for fb = 1:obj.Nmf
                    for pt = 1:obj.Npt
                        jittered_sfaps(:,pt,fb) = hr_shift_template(obj.sfaps(:,pt,fb), delays(fb)/obj.dt);
                    end
                end
                obj.muap = sum(jittered_sfaps, 3);
            else
                obj.muap = sum(obj.sfaps, 3);
            end
            
            muap = obj.muap;
        end
        
        %% Branching functions
        %function sim_nmj_branches_two_layers(obj, max_fibers_per_branch, endplate_area_center, branch_point_z_std, junction_z_std)
        function sim_nmj_branches_two_layers(obj, n_branches, endplate_area_center, branch_point_z_std, junction_z_std)
            
            % Cluster mfs to obtain branching point coordinates;
            obj.nerve_paths = zeros(obj.Nmf, 2);
            [idx, c] = kmeans(obj.mf_centers, n_branches);
            
            % Simulate axon terminal zs, this is z coordinate of branching point
            %std_from_area = mean(sqrt(diag(cov(obj.mf_centers))));
            obj.branch_points_xy = c;
            obj.branch_points_z = endplate_area_center + branch_point_z_std * randn(n_branches,1);
            
            obj.nmj_z = zeros(obj.Nmf,1);
            for i = 1:obj.Nmf
                %Simulate the junction z coordinate based on the branch coordinates
                obj.nmj_z(i) = obj.branch_points_z(idx(i)) + junction_z_std * randn();
            end
            
            obj.actual_center = [mean(obj.mf_centers), mean(obj.nmj_z)];
            for i = 1:obj.Nmf
                cluster_center = [obj.branch_points_xy(idx(i),:), obj.branch_points_z(idx(i))];
                nmj_coordinates = [obj.mf_centers(i,:), obj.nmj_z(i)];
                
                % Path from main nerve to cluster (or to first branching)
                obj.nerve_paths(i,1) = norm(obj.actual_center - cluster_center);
                % Path from first branching to the NMJ
                obj.nerve_paths(i,2) = norm(nmj_coordinates - cluster_center);
            end
            
        end
        
        
        
        function sim_nmj_branches_gaussian(obj, endplate_area_center, junction_z_std)
            % Cluster mfs to obtain branching point coordinates;
            obj.nerve_paths = zeros(obj.Nmf, 2);
            
            % Simulate axon terminal zs, this is z coordinate of branching point
            %std_from_area = mean(sqrt(diag(cov(obj.mf_centers)))); 
            
            obj.nmj_z = zeros(obj.Nmf,1);
            for i = 1:obj.Nmf
                %Simulate the junction z coordinate based on the branch coordinates
                obj.nmj_z(i) = endplate_area_center + junction_z_std * randn();
            end
            
            obj.actual_center = [mean(obj.mf_centers), mean(obj.nmj_z)];
            for i = 1:obj.Nmf
                cluster_center = obj.actual_center;
                nmj_coordinates = [obj.mf_centers(i,:), obj.nmj_z(i)];
                
                % Path from main nerve to cluster (or to first branching)
                obj.nerve_paths(i,1) = norm(nmj_coordinates - cluster_center);
                % Path from first branching to the NMJ
                obj.nerve_paths(i,2) = 0;
            end
            
        end            
           
        %% Spike train functions
        
        
        %% Force functions
        
        
        %% Other functions
        function distributed = distributeParameter(obj, par)
            if numel(par) == 1
                distributed = par*ones(obj.Nmf,1);
            else
                distributed = par;
            end
        end
        
    end
end


%         function sim_nmj_branches(obj)
%             obj.actual_center = mean(obj.mf_centers);
%             obj.nmj_b = zeros(obj.Nmf,1);
%             for i = 1:obj.Nmf
%                 obj.nmj_b(i) = sqrt(sum(([obj.mf_centers(i,:), obj.nmj_z(i)] - [obj.actual_center, obj.nmj_center]).^2));
%             end
%         end


%%
%            function sim_nmj_zs(obj, center, coeff_of_variation, dz)
%                 obj.nmj_center = center;
%
%                 if nargin < 3
%                     coeff_of_variation = 0.1;
%                 end
%
%                 obj.actual_center = mean(obj.mf_centers);
%                 std_from_area = mean(sqrt(diag(cov(obj.mf_centers))));
%
%                 %             if coeff_of_variation == 0
%                 %                 nmj_std = 0.5; %mm
%                 %             else
%                 nmj_std = std_from_area * coeff_of_variation;
%                 %             end
%
%                 if nargin == 4
%                     obj.nmj_z = dz * round((obj.nmj_center + nmj_std * randn(obj.Nmf,1))/dz);
%                 else
%                     obj.nmj_z = obj.nmj_center + nmj_std * randn(obj.Nmf,1);
%                 end
%            end

