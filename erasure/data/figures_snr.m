
par = 0.4:0.01:1;

data = load('snr_m.mat');
snr_m = transpose(data.snr_m);
data = load('snr_std.mat');
snr_std = transpose(data.snr_std);

hold on;

y = [snr_m(1:end-1) - snr_std(1:end-1); snr_m(2:end) - snr_std(2:end); flipud(snr_m(2:end) + snr_std(2:end));flipud(snr_m(1:end-1) + snr_std(1:end-1))];
x = [par(1:end-1); par(2:end);par(2:end);par(1:end-1)];

fill(x, y, 'blue','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');

errorbar(par(1:2:end), snr_m(1:2:end), snr_std(1:2:end),'-','LineWidth',1.2, 'CapSize',4, 'DisplayName','\Delta =10');


data = load('snr_m2.mat');
snr_m = transpose(data.snr_m);
data = load('snr_std2.mat');
snr_std = transpose(data.snr_std);

y = [snr_m(1:end-1) - snr_std(1:end-1); snr_m(2:end) - snr_std(2:end); flipud(snr_m(2:end) + snr_std(2:end));flipud(snr_m(1:end-1) + snr_std(1:end-1))];
x = [par(1:end-1); par(2:end);par(2:end);par(1:end-1)];

fill(x, y, 'red','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');

errorbar(par(1:2:end), snr_m(1:2:end), snr_std(1:2:end),'-','LineWidth',1.2, 'CapSize',4, 'DisplayName','\Delta =5');


data = load('snr_m3.mat');
snr_m = transpose(data.snr_m);
data = load('snr_std3.mat');
snr_std = transpose(data.snr_std);

y = [snr_m(1:end-1) - snr_std(1:end-1); snr_m(2:end) - snr_std(2:end); flipud(snr_m(2:end) + snr_std(2:end));flipud(snr_m(1:end-1) + snr_std(1:end-1))];
x = [par(1:end-1); par(2:end);par(2:end);par(1:end-1)];

fill(x, y, 'yellow','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');

errorbar(par(1:2:end), snr_m(1:2:end), snr_std(1:2:end),'-','LineWidth',1.2, 'CapSize',4, 'DisplayName','\Delta =2');

xlabel('$\mu$', 'Interpreter', 'latex');
ylabel('$\mathrm{SNR}$', 'Interpreter', 'latex');
legend('fontname','times');

