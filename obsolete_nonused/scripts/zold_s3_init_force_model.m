%% Twitch parameters 
%(see Fuglevand et al. - Models of Recruitment and Rate coding organisation in motor units, 1993.)
% Peak twitch forces: proportional to sizes
% Fuglevand's model:
%           twitch(t) = Pf/Tf * t .* exp(1 - t/Tf);
% Tf is time to peak;
% Pf is peak height;

% Peak heights
Pf = exp((log(rr)/N) * (1:N));

% Times to peak:
Tmax = 90/1000 * fs; % Maximum time-to-peak delay: 90 ms (Fuglevand)
Tr = 3; % Time to peak range
Tcoeff = log(rr)/log(Tr);

Tf = Tmax * (1./Pf).^(1/Tcoeff); % see Eq. 15 of Fuglevand

figure; plot(Pf); hold on; ylabel('P'); yyaxis right; plot(Tf); ylabel('T'); xlabel('MU');
title('P and T parameters distribution across the pool');

clear Tcoeff Tmax Tr

%% Generate twitches
max_twitch_len = ceil(5*max(Tf));
twitch_timeline = (0:(max_twitch_len-1))';
twitches_cell = cell(N,1);
twitches_mat = zeros(max_twitch_len, N);

for i = 1:N
    twitches_cell{i} = twitch_fuglevand(twitch_timeline, Pf(i), Tf(i));
    cut_tail = min(numel(twitches_cell{i}), ceil(5*Tf(i)));
    twitches_cell{i} = twitches_cell{i}(1:cut_tail);
    
    twitches_mat(:,i) = twitch_fuglevand(twitch_timeline, Pf(i), Tf(i));
end

timeline = ( 1:size(twitches_mat,1) )/fs;
figure; plot(timeline, twitches_mat, 'k');
title("Twitch waveforms"); xlabel("Time, s"); ylabel("Amplitude");

clear twitch_timeline cut_tail max_twitch_len



