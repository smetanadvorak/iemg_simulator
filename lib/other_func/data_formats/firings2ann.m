function [ filename ] = firings2ann( firings, filename, fs )
% [ ann ] = firings2ann( firings, filename, fs ) Writes firings to .ann file. If
% firings are sampled, specify 'fs' (sampling frequency) as well. 

if nargin < 2 
    filename = 'annotation';
    fs = nan;
end
if nargin < 3
    fs = nan;
end


if iscell(firings{1}) % if firings has multichannel structure (elements are cells)
    mkdir(filename);
    curdir = pwd; cd(filename);
    for i = 1:numel(firings)
        firings2ann(firings{i}, [filename, '_ch', num2str(i)], fs);
    end 
    cd(curdir);
    return
end



times = [];
mus = [];

for i = 1:length(firings)
    for j = 1:length(firings{i})
        times = [times; firings{i}(j)];
        mus = [mus; i];
    end
end

[times, ord] = sort(times);
if ~isnan(fs)
    times = times/fs;
end

mus = mus(ord);


dlmwrite([filename, '.ann'], [times mus], 'delimiter', ' ');





