function firings = shift_firings(firings, shifts, fs)
%function firings = shift_firings(firings, shifts)

if nargin < 3
    fs = 1;
end

if length(shifts) == 1
    shifts = shifts * ones(numel(firings),1);
end

if numel(firings) ~= numel(shifts)
    error('Number of motor neurons should be equal to length of shift vector');
end

for i = 1:numel(firings)
    firings{i} = firings{i} + shifts(i)/fs;
    firings{i}(firings{i}<=0) = [];
end


end

