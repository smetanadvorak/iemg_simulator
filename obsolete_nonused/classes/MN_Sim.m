classdef MN_Sim < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Main properties
        rt
        sz
        
        % Firing properties
        minfr
        maxfr
        frslope
        ipi_cv
        
        % Innervation properties
        area
        
        % MU
        MsU %Muscle unit (MU_Sim class)
        FU  %Force unit (to be implemented)
    end
    
    methods
        function obj = MN_Sim(sz, rt, minfr, maxfr, frslope, ipi_cv)
    
        end
        
        function [spikes, fr] = generate_spikes(excitation)
            
            T = length(excitation);
            spikes = zeros(T,1);
            fr = zeror(T,1);
            ready_to_fire = 0;
            for t = 1:T
                fr(t) = calculate_fr(obj.rt(m), obj.minfr(m), obj.maxfr(m), obj.frs(m), excitation(t));
                if excitation(t) > obj.rt
                    %Plan next firing
                    if ~ready_to_fire
                        ipi = fs/fr(t);
                        ipi = ipi + randn * ipi / obj.ipi_cv;
                        next_firing = t + round(ipi);
                        ready_to_fire = 1;
                    elseif t == next_firing
                        spikes(t) = 1;
                        ready_to_fire = 0;
                    end
                    
                    if excitation(t) < obj.rt(m)
                        ready_to_fire = 0;
                        next_firing = nan;
                    end
                end
            end
            
        end
        
        function prob = calc_adoption_probability()
        
        
        end
    end
    
end

