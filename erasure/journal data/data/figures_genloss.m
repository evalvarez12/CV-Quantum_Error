clear all;
close all;
set(0,'defaultTextInterpreter','latex');

% PD
N = 2e5;
Nsubset = 2e3;

par = linspace(0.4, 1, 60);
sigma_coh = 10;
F_class = (1 + (1/sigma_coh))/(2 + (1/sigma_coh));

X = repmat(par, Nsubset, 1);
X = transpose(X);

data = load('F_noerased_V3_sigma10.mat');
Fm = data.Fm;
Fmean = mean(transpose(Fm));
Ncl = sum(transpose(Fm) > F_class)/N;

% data = load('F_pd_std.mat');
% Fstd = data.Fstd;

data = load('F_dir_sigma10.mat');
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

yyaxis left
% title('Plots with Different y-Scales')
xlabel('$T$', 'Interpreter', 'latex');


plot(par, Fmean, 'LineWidth', 1.7, 'DisplayName', 'Protocol');
plot(par, Fmean_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');

yyaxis right
ylabel('$\%> F_\mathrm{class}$', 'Interpreter', 'latex');

plot(par, Ncl, 'LineWidth', 1.7, 'DisplayName', 'Protocol');
plot(par, Ncl_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');
legend('fontname','times', 'Location','southeast');

subplot(2, 2, 2)
title('Protocol')  
% figure;
hold on;
scatter(X, Fm(:, 1:Nsubset), 1, 'MarkerFaceAlpha',.01,'MarkerEdgeAlpha',.01);
plot(par, Fmean, 'LineWidth', 1.7);
% plot(x, Ncl, 'LineWidth', 1.7);
xlabel('$\sigma$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');

plot(par, class, 'LineWidth', 1, 'Color', 'red');

subplot(2, 2,4)
title('Direct')  
% figure;
hold on;
scatter(X, Fm_dir(:, 1:Nsubset), 1, '+' ,'MarkerFaceAlpha',.01,'MarkerEdgeAlpha',.01);
plot(par, Fmean_dir, 'LineWidth', 1.7);
plot(par, class, 'LineWidth', 1, 'Color', 'red');

% plot(x, Ncl_dir, 'LineWidth', 1.7);
xlabel('$\sigma$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');

