function h = get_elementary_current_response( z, z_electrode, r, sigma_r, sigma_z )

    if nargin < 5
        sigma_z = 0.33 * 1000; %S/mm
        %values from Andreassen and Rosenfalck 1980
    end
    if nargin < 4
        sigma_r = 0.063 * 1000; %S/mm
        %values from Andreassen and Rosenfalck 1980
    end
        
    h = 1/4/pi/sigma_r ./ sqrt(sigma_z/sigma_r * r.^2 + (z-z_electrode).^2); % SIMULATION OF MACRO EMG MOTOR UNIT POTENTIALS, SANJEEV NANDEDKAR and ERIK STALBERG
    %h = 1 ./ sqrt((sigma_r^2) * r^2 +  (sigma_z^2) * (z-z_electrode).^2); % My logic
end


