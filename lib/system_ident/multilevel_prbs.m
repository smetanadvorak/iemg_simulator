function sig = multilevel_prbs(T, period, n_levels)

levels = linspace(0, 1, n_levels);
warning('off','Ident:dataprocess:idinput7');
sig = idinput(2*T*numel(levels), 'prbs', [0 1/period], [0,1]);
warning('on','Ident:dataprocess:idinput7');

sig = sum(reshape(sig, [], numel(levels)),2);
sig = majority_filter(sig, 2*round(period/8)+1)/numel(levels);


end

