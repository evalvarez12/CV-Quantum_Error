clear all;
close all;
set(0,'defaultTextInterpreter','latex');

% newcolors = [0.9290, 0.6940, 0.1250
%              0.8500, 0.3250, 0.0980
%              0, 0.4470, 0.7410
%              0.25 0.80 0.54];
%          
% colororder(newcolors)


% PD
N = 10000;

rec = '0.2';
dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000];

% rec = '0.15';
% dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600];

% rec = '0.1';
% dists = [1000 1200 1400 1600 1800 2000 2200 2400];

sigma_coh = 10;
F_class = (1 + (1/sigma_coh))/(2 + (1/sigma_coh));


data = load(['Fscatter_erasure_phasesc_rec=', rec,'_V3.mat']);
Fm = data.Fm;
Fmean = mean(transpose(Fm));
Fstd = std(transpose(Fm));

data2 = load(['Fscatter_erasure_phasesc_rec=', rec, '_V10.mat']);
Fm2 = data2.Fm;
Fmean2 = mean(transpose(Fm2));
Fstd2 = std(transpose(Fm2));


% data3 = load(['Fscatter_erasure_phasesc_rec=', rec, '_V6.mat']);
% Fm3 = data3.Fm;
% Fmean3 = mean(transpose(Fm3));
% Fstd3 = std(transpose(Fm3));


% data = load('F_pd_std.mat');
% Fstd = data.Fstd;

data = load(['Fscatter_phasesc_rec=', rec, '_V10.mat']);
Fm_dir = data.Fm;
Fmean_dir = mean(transpose(Fm_dir));
Fstd_dir = std(transpose(Fm_dir));

y = [Fmean(1:end-1) - Fstd(1:end-1); Fmean(2:end) - Fstd(2:end); flipud(Fmean(2:end) + Fstd(2:end)); flipud(Fmean(1:end-1) + Fstd(1:end-1))];

y2 = [Fmean2(1:end-1) - Fstd2(1:end-1); Fmean2(2:end) - Fstd2(2:end); flipud(Fmean2(2:end) + Fstd2(2:end)); flipud(Fmean2(1:end-1) + Fstd2(1:end-1))];

x = [dists(1:end-1); dists(2:end); dists(2:end); dists(1:end-1)];

y_dir = [Fmean_dir(1:end-1) - Fstd_dir(1:end-1); Fmean_dir(2:end) - Fstd_dir(2:end); flipud(Fmean_dir(2:end) + Fstd_dir(2:end)); flipud(Fmean_dir(1:end-1) + Fstd_dir(1:end-1))];

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
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');


fill(x, y_dir, 'black','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');

fill(x, y, 'blue','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');

fill(x, y2, 'red','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');

plot(dists, Fmean, 's-', 'LineWidth', 1.2, 'DisplayName', 'Erasure - V=3');

plot(dists, Fmean2, 'v-', 'LineWidth', 1.2, 'DisplayName', 'Erasure - V=10');

% plot(par(1:10:end), Fmean2(1:10:end),'*-', 'LineWidth', 1.2, 'DisplayName', 'V=6');

plot(dists, Fmean_dir, 'k-', 'LineWidth', 1.7, 'DisplayName', 'No erasure');


% xlim([par(1) par(end)])
ylim([.5 1])
legend('Location','northeast')


savefigures('scatter_phasesc_erasure_02');



% yyaxis right
% ylabel('$\%> F_\mathrm{class}$', 'Interpreter', 'latex');
% 
% plot(par, Ncl, 'LineWidth', 1.7, 'DisplayName', 'Protocol');
% plot(par, Ncl_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');
% legend('fontname','times', 'Location','southeast');


