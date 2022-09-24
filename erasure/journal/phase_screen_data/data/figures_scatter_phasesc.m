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
N = 200;
rec = '0.1';

dists = [1000 1400 1800 2200 2600 3000];
sigma_coh = 10;
F_class = (1 + (1/sigma_coh))/(2 + (1/sigma_coh));


data = load(['Fscatter_phasesc_rec=', rec,'_V3.mat']);
Fm = data.Fm;
Fmean = mean(transpose(Fm));


data2 = load(['Fscatter_phasesc_rec=', rec, '_V10.mat']);
Fm2 = data2.Fm;
Fmean2 = mean(transpose(Fm2));




% data = load('F_pd_std.mat');
% Fstd = data.Fstd;

data = load(['Fscatter_phasesc_rec=', rec, '_dir.mat']);
Fm_dir = data.Fm_dir;
Fmean_dir = mean(transpose(Fm_dir));


% data = load('F_pd_std_dir.mat');
% Fstd_dir = data.Fstd_dir;

class = ones(size(dists))*F_class;

set(gca,'fontname','times') 
% subplot(2, 2, [1 3])
% figure;
hold on;

% yyaxis left
% title('Plots with Different y-Scales')
xlabel('$T$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');

plot(dists, Fmean2, 'v-', 'LineWidth', 1.2, 'DisplayName', 'V=3');

plot(dists, Fmean, '*-', 'LineWidth', 1.2, 'DisplayName', 'V=10');

% plot(par(1:10:end), Fmean2(1:10:end),'*-', 'LineWidth', 1.2, 'DisplayName', 'V=6');

plot(dists, Fmean_dir, 'k-', 'LineWidth', 1.7, 'DisplayName', 'Direct');


% xlim([par(1) par(end)])
% ylim([.55 .95])
legend('Location','northwest')





% yyaxis right
% ylabel('$\%> F_\mathrm{class}$', 'Interpreter', 'latex');
% 
% plot(par, Ncl, 'LineWidth', 1.7, 'DisplayName', 'Protocol');
% plot(par, Ncl_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');
% legend('fontname','times', 'Location','southeast');


