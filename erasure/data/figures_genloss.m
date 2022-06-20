clear all;
close all;

% PD
N = 2e5;
Nsubset = 2e3;

% mu
x = 0.4:0.002:1;
% sigma
% x = 0:0.001:.3;

sigma_coh = 10;
F_class = (1 + (1/sigma_coh))/(2 + (1/sigma_coh));

X = repmat(x, Nsubset, 1);
X = transpose(X);

data = load('F_mu.mat');
Fm = data.Fm;
Fmean = mean(transpose(Fm));
Ncl = sum(transpose(Fm) > F_class)/N;

% data = load('F_pd_std.mat');
% Fstd = data.Fstd;

data = load('F_mu_dir.mat');
Fm_dir = data.Fm_dir;
Fmean_dir = mean(transpose(Fm_dir));
Ncl_dir = sum(transpose(Fm_dir) > 0.5)/N;

% data = load('F_pd_std_dir.mat');
% Fstd_dir = data.Fstd_dir;

class = ones(size(x))*F_class;

set(gca,'fontname','times') 
subplot(2, 2, [1 3])
% figure;
hold on;

yyaxis left
% title('Plots with Different y-Scales')
xlabel('$\mu$', 'Interpreter', 'latex');


plot(x, Fmean, 'LineWidth', 1.7, 'DisplayName', 'Code');
plot(x, Fmean_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');
ylabel('$\mathcal{F}$', 'Interpreter', 'latex');

yyaxis right
ylabel('$\%$', 'Interpreter', 'latex');

plot(x, Ncl, 'LineWidth', 1.7, 'DisplayName', 'Code');
plot(x, Ncl_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');
legend('fontname','times', 'Location','southwest');

subplot(2, 2, 2)
title('Code')  
% figure;
hold on;
scatter(X, Fm(:, 1:Nsubset), 1, 'MarkerFaceAlpha',.01,'MarkerEdgeAlpha',.01);
plot(x, Fmean, 'LineWidth', 1.7);
% plot(x, Ncl, 'LineWidth', 1.7);
xlabel('$\mu$', 'Interpreter', 'latex');
ylabel('$\mathcal{F}$', 'Interpreter', 'latex');

plot(x, class, 'LineWidth', 1, 'Color', 'red');

subplot(2, 2,4)
title('Direct')  
% figure;
hold on;
scatter(X, Fm_dir(:, 1:Nsubset), 1, '+' ,'MarkerFaceAlpha',.01,'MarkerEdgeAlpha',.01);
plot(x, Fmean_dir, 'LineWidth', 1.7);
% plot(x, Ncl_dir, 'LineWidth', 1.7);
xlabel('$\mu$', 'Interpreter', 'latex');
ylabel('$\mathcal{F}$', 'Interpreter', 'latex');

plot(x, class, 'LineWidth', 1, 'Color', 'red');




% errorbar(X, Fm, Fstd);
% errorbar(X, Fm_dir, Fstd_dir);

% % T0
% X = 0:0.02:0.8;
% 
% data = load('F_T0_mean.mat');
% Fm = data.Fm;
% 
% data = load('F_T0_std.mat');
% Fstd = data.Fstd;
% 
% data = load('F_T0_mean_dir.mat');
% Fm_dir = data.Fm_dir;
% 
% data = load('F_T0_std_dir.mat');
% Fstd_dir = data.Fstd_dir;
% 
% figure;
% hold on;
% errorbar(X, Fm, Fstd);
% errorbar(X, Fm_dir, Fstd_dir);
% 
% % sigma
% X = 0:0.002:0.2;
% 
% data = load('F_sigma_mean.mat');
% Fm = data.Fm;
% 
% data = load('F_sigma_std.mat');
% Fstd = data.Fstd;
% 
% data = load('F_sigma_mean_dir.mat');
% Fm_dir = data.Fm_dir;
% 
% data = load('F_sigma_std_dir.mat');
% Fstd_dir = data.Fstd_dir;
% 
% figure;
% hold on;
% errorbar(X, Fm, Fstd);
% errorbar(X, Fm_dir, Fstd_dir);