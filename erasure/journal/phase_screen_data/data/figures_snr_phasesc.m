clear all;
close all;
set(0,'defaultTextInterpreter','latex');



rec = '0.2';
dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000];

% rec = '0.15';
% dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600];

% rec = '0.1';
% dists = [1000 1200 1400 1600 1800 2000 2200 2400];

sigma_coh = 10;
F_class = (1 + (1/sigma_coh))/(2 + (1/sigma_coh));


data = load(['snrm_phasesc_rec=', rec,'_delta=5.mat']);
snr_m = transpose(data.snr_m);
data = load(['snrstd_phasesc_rec=', rec,'_delta=5.mat']);
snr_std = transpose(data.snr_std);


data = load(['snrm_phasesc_rec=', rec,'_delta=10.mat']);
snr_m2 = transpose(data.snr_m);
data = load(['snrstd_phasesc_rec=', rec,'_delta=10.mat']);
snr_std2 = transpose(data.snr_std);





% data = load('F_pd_std.mat');
% Fstd = data.Fstd;

% data = load(['Fscatter_phasesc_rec=', rec, '_dir.mat']);
% Fm_dir = data.Fm_dir;
% Fmean_dir = mean(transpose(Fm_dir));
% Fstd_dir = std(transpose(Fm_dir));

y = [snr_m(1:end-1) - snr_std(1:end-1); snr_m(2:end) - snr_std(2:end); flipud(snr_m(2:end) + snr_std(2:end)); flipud(snr_m(1:end-1) + snr_std(1:end-1))];

y2 = [snr_m2(1:end-1) - snr_std2(1:end-1); snr_m2(2:end) - snr_std2(2:end); flipud(snr_m2(2:end) + snr_std2(2:end)); flipud(snr_m2(1:end-1) + snr_std2(1:end-1))];

x = [dists(1:end-1); dists(2:end); dists(2:end); dists(1:end-1)];

% y_dir = [Fmean_dir(1:end-1) - Fstd_dir(1:end-1); Fmean_dir(2:end) - Fstd_dir(2:end); flipud(Fmean_dir(2:end) + Fstd_dir(2:end)); flipud(Fmean_dir(1:end-1) + Fstd_dir(1:end-1))];

% data = load('F_pd_std_dir.mat');
% Fstd_dir = data.Fstd_dir;

class = ones(size(dists))*F_class;

set(gca,'fontname','times') 
% subplot(2, 2, [1 3])
% figure;
hold on;

% yyaxis left
% title('Plots with Different y-Scales')
xlabel('$L[m]$', 'Interpreter', 'latex');
ylabel('$\mathrm{SNR}$', 'Interpreter', 'latex');

plot(dists, snr_m, 'v-', 'LineWidth', 1.2, 'DisplayName', '\Delta=5');

plot(dists, snr_m2, '*-', 'LineWidth', 1.2, 'DisplayName', '\Delta=10');



fill(x, y, 'blue','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');

fill(x, y2, 'red','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');

% xlim([par(1) par(end)])
% ylim([.5 1])
legend('Location','northeast')


% savefigures('scatter_phasesc_02');



% yyaxis right
% ylabel('$\%> F_\mathrm{class}$', 'Interpreter', 'latex');
% 
% plot(par, Ncl, 'LineWidth', 1.7, 'DisplayName', 'Protocol');
% plot(par, Ncl_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');
% legend('fontname','times', 'Location','southeast');


