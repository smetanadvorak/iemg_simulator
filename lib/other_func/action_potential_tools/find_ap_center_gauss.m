function [centers, widths] = find_ap_center_gauss( ap, winlen )
% Find the center of an action potential using gaussian filter over
% absolute values
len = size(ap, 1);
nch = size(ap, 2);

if nargin < 2
    winlen = len;
end

gw = gausswin(winlen, 3);
centers = zeros(nch,1);
widths = zeros(nch, 1);

for ch = 1:nch
    filtered = conv(abs(ap(:,ch)), gw, 'same');
    [pks, locs, wid] = findpeaks(filtered);
    
    if numel(pks) > 1
        warning('%d possible centers obtained during action potential centering, closest to the middle is chosen', numel(pks));
        [~, closest] = min(locs - len/2);
        pks = pks(closest);
        locs = locs(closest);
        wid = wid(closest);
    end
    
    centers(ch) = locs;
    widths(ch) = wid;
end
end

