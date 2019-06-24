function cv = assign_mf_cv( diameters, max_cv, intercept, min_diameter, max_diameter )

if nargin < 5
    max_diameter = 85e-6; %um, %Taken from S. D. Nandedkar E. Stalberg
                            %Simulation of single muscle fibre action potentials by
end

if nargin < 4
    min_diameter = 22e-6; %um, % From the same source
end

if nargin < 3
    intercept = 2200; %mm/s 	% From the same source
end

if nargin < 2
    max_cv = 5200; %mm/s      % From the same source
end

diameters(diameters < min_diameter) = min_diameter;
diameters(diameters > max_diameter) = max_diameter;


slope = (max_cv-intercept)/(max_diameter-min_diameter); % m/s/um 

cv = intercept + (diameters-min_diameter) * slope;

end

