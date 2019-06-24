function [sig_filt] = majority_filter(class_flow, win)
% Do a majority voting on discrete-value sequence (kind of a erosion),
% values should be integers

if ~mod(win,2)
    error('Window should be of odd length');
end

%sig = round(sig);

base = min(class_flow) + 1;
class_flow = class_flow + base;

vals = zeros(numel(unique(class_flow)), 1);

sig_filt = class_flow;
hp = ceil(win/2);
hm = floor(win/2);

for i = 1:win
    vals(class_flow(i)) = vals(class_flow(i)) + 1;
end

for i = hm+1 : length(class_flow)-hp
    vals(class_flow(i-hm)) = vals(class_flow(i-hm)) - 1;
    vals(class_flow(i+hp)) = vals(class_flow(i+hp)) + 1;
    [~, max_vote] = max(vals);
    sig_filt(i) = max_vote;
end

sig_filt = sig_filt - base;
end

