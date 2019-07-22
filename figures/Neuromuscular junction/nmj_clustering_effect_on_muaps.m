to_take_time = 15:250;
timeline = (to_take_time - to_take_time(1))/fs*1000; %ms

analysisN = 84;
figure; 


for i = 1:size(MUs(analysisN).sfaps, 3)
    plot(timeline, squeeze(MUs(analysisN).sfaps(to_take_time,:,i)) * electrode_diff_mat(1,:)','k-'); hold on;
end
ylabel('SFAPs');

yyaxis right
plot(timeline, MUs(analysisN).muap(to_take_time,:) * electrode_diff_mat(1,:)', 'linewidth', 2);
%plot(sum(MUs(mu).sfaps, 3) * electrode_diff_mat(1,:)', '--');
ylabel('MUAP');

%%align zero for left and right
yyaxis right; ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
yyaxis left; yliml = get(gca,'Ylim');
if yliml(2)*ratio<yliml(1)
    set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
else
    set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
end

xlim([timeline(1), timeline(end)]);

xlabel('Time, ms');
saveas(gcf, 'MUAPs stucture', 'png');