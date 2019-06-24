classdef MF_Sim < handle
    
    properties
        % Important constants
        
        
        % Geometry
        limits; % 2x1 vector with left and right fiber ends (z axis)
        center; % coordinates of the mf center (x,y)
        nmj_z;
        
        % Electrical characteristics
        diameter;
        cv;
        
        % Motor unit assignment
        mu;
        mn;
        
        
        
    end
    
    methods
        function obj = MF_sim(limits, center, diameter)
            obj.limits = limits;
            obj.center = center;
            obj.diameter = diameter;
        end
        
        function assign_cv(obj)
            obj.cv = assign_mf_cv(obj.diameter);
        end
        
        function current_density = calc_current_density(obj, dt, dz)
            t = 0:dt:(2 * max([(obj.nmj_z - obj.limits(1))./obj.cv ; (obj.limits(2) - obj.nmj_z)./obj.cv]));
            t = t';
 
            % Symmetrical z relative to the nmj coordinate (gives symmetrical spatial current densities):
            z_left = obj.nmj_z(fb):-dz:(obj.limits(1));
            z_right = obj.nmj_z(fb):dz:(obj.limits(2));
            z = transpose([z_left(end:-1:1), z_right(2:end)]);
            
            % Get current density
            current_density = get_current_density_symm(t, z, obj.nmj_z, obj.limits(1) - obj.nmj_z, obj.nmj_z - obj.limits(2), obj.cv, obj.diameter);
        end
        
        function electric_distance = calc_electric_distance(obj, pt)
            electric_distance = get_electric_distance(z, pt(3), sqrt(sum((pt(1:2) - obj.center).^2)));
        end
        
        function cd = get_current_density_fast(obj, precalculated)
        end
          
    end
    
end

