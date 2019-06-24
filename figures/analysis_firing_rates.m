%% Calculate firing rates
FR = zeros(T,N);
for m = 1:N
    for i = 1:T
        FR(i,m) = calculate_fr(rt(m), minfr(m), maxfr(m), frs(m), excitation(i));
    end
end

figure; plot(emg_timeline, FR); 
ylim([min(minfr), max(maxfr)]);
ylabel('Firing rates, pulses/s');
yyaxis right;
plot(emg_timeline, excitation, 'linewidth', 2);
ylabel('Excitation, normalized units');
xlabel('Time, sec');
title('Motor neurons'' firining rates and excitation vs time');