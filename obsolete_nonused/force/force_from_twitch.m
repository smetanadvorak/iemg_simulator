function [ force ] = force_from_twitch( Tf, Pf, sawtooth, win_len )

force = zeros(size(sawtooth,1),1);

for i = 1:size(force,1)    
    delays = sawtooth(i,:)';
    force(i) = sum(twitch_texpt(delays(:), Tf(:), Pf(:), win_len));
end

end

