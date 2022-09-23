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
N = 5e5;
Nsubset = 5e2;

par = .5:0.002:1;
sigma_coh = 10;
F_class = (1 + (1/sigma_coh))/(2 + (1/sigma_coh));

X = repmat(par, Nsubset, 1);
X = transpose(X);


sigma = '05';

data = load(['Fscatter_V10_sigma', sigma, '.mat']);
Fm = data.Fm;
Fmean = mean(transpose(Fm));
Ncl = sum(transpose(Fm) > F_class)/N;


data2 = load(['Fscatter_V10_sigma', sigma, '.mat']);
Fm2 = data2.Fm;
Fmean2 = mean(transpose(Fm2));

data3 = load(['Fscatter_V3_sigma', sigma, '.mat']);
Fm3 = data3.Fm;
Fmean3 = mean(transpose(Fm3));


% data = load('F_pd_std.mat');
% Fstd = data.Fstd;

data = load(['Fscatter_dir_sigma', sigma, '.mat']);
Fm_dir = data.Fm_dir;
Fmean_dir = mean(transpose(Fm_dir));
Ncl_dir = sum(transpose(Fm_dir) > 0.5)/N;

% data = load('F_pd_std_dir.mat');
% Fstd_dir = data.Fstd_dir;

class = ones(size(par))*F_class;

set(gca,'fontname','times') 
subplot(2, 2, [1 3])
% figure;
hold on;

% yyaxis left
% title('Plots with Different y-Scales')
xlabel('$T$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');

plot(par(1:10:end), Fmean3(1:10:end), 'v-', 'LineWidth', 1.2, 'DisplayName', 'V=3');

plot(par(1:10:end), Fmean(1:10:end), '*-', 'LineWidth', 1.2, 'DisplayName', 'V=10');

% plot(par(1:10:end), Fmean2(1:10:end),'*-', 'LineWidth', 1.2, 'DisplayName', 'V=6');

plot(par(1:10:end), Fmean_dir(1:10:end), 'k-', 'LineWidth', 1.7, 'DisplayName', 'Direct');


xlim([par(1) par(end)])
ylim([.55 .95])
legend('Location','northwest')





% yyaxis right
% ylabel('$\%> F_\mathrm{class}$', 'Interpreter', 'latex');
% 
% plot(par, Ncl, 'LineWidth', 1.7, 'DisplayName', 'Protocol');
% plot(par, Ncl_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');
% legend('fontname','times', 'Location','southeast');





subplot(2, 2, 2)
title('Protocol')  
% figure;
hold on;
scatter(X, Fm(:, 1:Nsubset), 1, 'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.01);
plot(par, Fmean, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.7);
% plot(x, Ncl, 'LineWidth', 1.7);
xlabel('$\mu$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');

plot(par, class, 'LineWidth', 1, 'Color', 'red');

xlim([par(1) par(end)])

subplot(2, 2, 4)
title('Direct')  
% figure;
hold on;
scatter(X, Fm_dir(:, 1:Nsubset), 1, '+' ,'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.01);
plot(par,  Fmean_dir, 'Color', [0, 0.4470, 0.7410],  'LineWidth', 1.7);
plot(par, class, 'LineWidth', 1, 'Color', 'red');

% plot(x, Ncl_dir, 'LineWidth', 1.7);
xlabel('$\mu$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');

xlim([par(1) par(end)])


% savefigures2('scatter.pdf')