classdef Electrode < handle
    properties
        pts_origin; %Facet centers in the electrode's frame
        normals_origin; %Facet orientation in the electrode's frame, (set to [0,0,0] for point model)
        
        pts_init;   %Facet centers in initial position in muscle frame
        normals_init;
        
        pts;        %Facet centers along trajectory in muscle frame
        normals;    %Facet orientation in muscle frame ([0,0,0] for point model)
        
        diff_mat;    %size(pts,1) x 1 difference vector for difference mode electrodes
        
        n_points;
        n_channels;
        
        pos;        %Electrode's origin position in muscle frame
        orn;        %Electrode's orientation in muscle frame (roll, pitch, yaw), [rad]
        
        type;
        type_short;
        
        simplified_pts; % For visualisation, if too many facets
        
        traj_transforms;
        traj_step;
        n_nodes;
        traj_mixing_fun;
        traj_mixing_mat;
        
    end
    
    methods
        function obj = Electrode()
            
            
        end
        
        function obj = set_position(obj, pos, orn)
            obj.orn = orn(:);
            if isempty(obj.orn)
                obj.orn = [0;0;0];
            end
            
            obj.pos = pos(:);
            
            % Rotate in muscle frame
            obj.pts_init = obj.pts_origin;
            obj.pts_init = rodrigues_rot(obj.pts_init, [1, 0, 0],  obj.orn(1));
            obj.pts_init = rodrigues_rot(obj.pts_init, [0, 1, 0],  obj.orn(2));
            obj.pts_init = rodrigues_rot(obj.pts_init, [0, 0, 1],  obj.orn(3));
            
            % Translate in muscle frame
            obj.pts_init = obj.pts_init + repmat(obj.pos', size(obj.pts_init,1), 1);
            obj.pts = obj.pts_init;
        end
        
        
        function set_linear_trajectory(obj, distance, n_nodes)
            % Sets linear trajectory along the electrode's axis
            
            if nargin < 3
                n_nodes = max(ceil(distance/0.5), 1); %mm
            end
            
            obj.n_nodes = n_nodes;
            obj.traj_step = distance/obj.n_nodes;
            
            obj.traj_transforms = linspace(0, distance, obj.n_nodes);
            obj.traj_transforms = [ zeros(length(obj.traj_transforms), 2), ...
                obj.traj_transforms(:), ...
                zeros(length(obj.traj_transforms), 3)];
            
            obj.traj_transforms(:,1:3) = rodrigues_rot(obj.traj_transforms(:,1:3), [1, 0, 0],  obj.orn(1));
            obj.traj_transforms(:,1:3) = rodrigues_rot(obj.traj_transforms(:,1:3), [0, 1, 0],  obj.orn(2));
            obj.traj_transforms(:,1:3) = rodrigues_rot(obj.traj_transforms(:,1:3), [0, 0, 1],  obj.orn(3));
            
            obj.calc_observation_points();
        end
        
        function calc_observation_points(obj)
            %% Calculate all observation points along the trajectory:
            obj.pts = [];
            for i = 1:obj.n_nodes
                snapshot = rotate_and_translate(obj.pts_init, obj.traj_transforms(i, 4:6), obj.traj_transforms(i, 1:3));
                obj.pts = [obj.pts; snapshot];
            end
            
            % Extended electrode matrix
            obj.diff_mat = repmat(obj.diff_mat, 1, obj.n_nodes);
            
            % Position-dependent mixing matrix to approximate the simulation output
            % between the nodes.
            % Parameter t is bounded between 0 (initial position) and 1 (finish of the
            % trajectory). You can map any behaviour into that parameter, for example,
            % the contraction force, or time.
            
            %% What makes the electrode translate: time or force.
            % Normalized parameter (between 0 and 1);
            obj.traj_mixing_fun = @(t, n_nodes, node) max(0, 1 - (n_nodes-1) * abs(t - (node-1)/(n_nodes-1+eps)) );
            obj.traj_mixing_mat = @(t, n_nodes, n_channels) diag(reshape(repmat(obj.traj_mixing_fun(t, n_nodes, 1:n_nodes), obj.n_points, 1), [], 1));
            
        end
        
        function draw(obj, ax, whatplot)
            if nargin < 2
                figure;
                ax = axes;
            end
            if nargin < 3
                whatplot = 'init';
            end
            
            switch whatplot
                case 'init'
                    plot3(ax, obj.pts_init(:,1),obj.pts_init(:,2),obj.pts_init(:,3),'*b');
                case 'trajectory'
                    plot3(ax, obj.pts);
            end
        end
    end
end


%         function set_trajectory(obj, nodes, orns)
%             warning('\n Method not implemented!\n');
%             return

%             % Electrode is modeled as a rigid point cloud, so that its own geometry does
%             % not change while moving (which can be wrong in case of fine-wire
%             % electrodes, so pay attention to that).
%             % The trajectory is described as a number o transformations from
%             % the initial node to the current one, it contains six numbers: x, y, z
%             % and alpha, beta, gamma (yaw, pitch and roll angles).
%
%             if nargin < 2 || isempty(nodes) || all(nodes == 0)
%                 % Stationary electrode:
%                 obj.traj_transforms =    [0,0,0,0,0,0];
%             else
%                 obj.traj_transforms = nodes;
%                 if nargin > 3
%                     obj.orns = orns;
%                 else
%                     obj.orns = repmat(obj.orn(:)', size(nodes, 1), 1);
%                 end
%             end
%
%             %% Slight shifting along x axis (transversally to the fibers)
%             % traj_step = 0.25;
%             % traj_transforms = linspace(0, 1, round(1/traj_step));
%             % traj_transforms = [traj_transforms(:), zeros(length(traj_transforms), 5)];
%
%             %% Slight shifting along z axis (longitudinally to the fibers)
%             % traj_step = 0.25;
%             % traj_transforms = linspace(0, 1, round(1/traj_step));
%             % traj_transforms = [zeros(length(traj_transforms), 2), traj_transforms(:), zeros(length(traj_transforms), 3)];
%
%             %% Scan across muscle
%             % traj_step = 0.25; %mm
%             % traj_transforms =  -Rmuscle : traj_step : Rmuscle; % Translate along x Axis
%             % traj_transforms = [traj_transforms(:), zeros(length(traj_transforms), 5)];
%
%             %% Direct setting
%             % traj_transforms =    [0.00,0,0,0,0,0;
%             %                       0.10,0,0,0,0,0;
%             %                       0.20,0,0,0,0,0;
%             %                       0.30,0,0,0,0,0;
%             %                       0.40,0,0,0,0,0;
%             %                       0.50,0,0,0,0,0]/2; %x,y,z, yaw,pitch,roll;
%
%             %% Other trajectories:

%         end

