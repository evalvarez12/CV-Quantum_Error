clear all;
close all;
set(0,'defaultTextInterpreter','latex');

par = linspace(0.4, 1, 24);

sigma_coh = 10;
F_class = (1 + (1/sigma_coh))/(2 + (1/sigma_coh));


% dataV1 = load('F_compare_V1_sigma10.mat');
dataV1 = load('F_compare_V1_sigma10_sigma0.mat');
% dataV1 = load('F_compare_V1_sigma10_sigma05.mat');

dataV1 = dataV1.results;

F1V1 = dataV1(:, :, 2);
Fm_1V1= mean(F1V1, 2);


dataV3 = load('F_compare_V3_sigma10_sigma10.mat');
% dataV3 = load('F_compare_V3_sigma10_sigma0.mat');
% dataV3 = load('F_compare_V3_sigma10_sigma05.mat');

dataV3 = dataV3.results;

FnoV3 = dataV3(:, :, 1);
F1V3 = dataV3(:, :, 2);
F2V3 = dataV3(:, :, 3);
F3V3 = dataV3(:, :, 4);
FdirV3= dataV3(:, :, 5);

Fm_noV3 = mean(FnoV3, 2);
Fm_1V3= mean(F1V3, 2);
Fm_2V3= mean(F2V3, 2);
Fm_3V3= mean(F3V3, 2);
Fm_dirV3 =  mean(FdirV3, 2);


dataV5 = load('F_compare_V5_sigma10_sigma10.mat');
% dataV5 = load('F_compare_V5_sigma10_sigma0.mat');
% dataV5 = load('F_compare_V5_sigma10_sigma05.mat');

dataV5 = dataV5.results;

F1V5 = dataV5(:, :, 2);
Fm_1V5= mean(F1V5, 2);

dataV10 = load('F_compare_V10_sigma10_sigma10.mat');
% dataV10 = load('F_compare_V10_sigma10_sigma0.mat');
% dataV10 = load('F_compare_V10_sigma10_sigma05.mat');

dataV10 = dataV10.results;

F1V10 = dataV10(:, :, 2);
Fm_1V10= mean(F1V10, 2);

dataV20 = load('F_compare_V20_sigma10_sigma10.mat');
% dataV20 = load('F_compare_V20_sigma10_sigma0.mat');
% dataV20 = load('F_compare_V20_sigma10_sigma05.mat');

dataV20 = dataV20.results;

F1V20 = dataV20(:, :, 2);
Fm_1V20= mean(F1V20, 2);

FnoV20 = dataV20(:, :, 1);
Fm_noV20 = mean(FnoV20, 2);

Fdir = dataV20(:, :, 5);
Fm_dir = mean(Fdir, 2);


set(gca,'fontname','times') 
hold on;

plot(par, Fm_noV3, 'LineWidth', 1.7, 'DisplayName', 'No erasure');
plot(par, Fm_1V3, 'LineWidth', 1.7, 'DisplayName', 'Erasure - $V=3$');
plot(par, Fm_1V5, 'LineWidth', 1.7, 'DisplayName', 'Erasure - $V=5$');
plot(par, Fm_1V10, 'LineWidth', 1.7, 'DisplayName', 'Erasure - $V=10$');
plot(par, Fm_1V20, 'LineWidth', 1.7, 'DisplayName', 'Erasure - $V=20$');
% plot(par, Fm_1V1, 'LineWidth', 1.7, 'DisplayName', 'Erasure - $V=1$');

% plot(par, Fm_noV20, 'LineWidth', 1.7, 'DisplayName', 'No erasure');
% plot(par, Fm_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');

% 
% plot(par, Fm_2, 'LineWidth', 1.7, 'DisplayName', '2 erased');
% plot(par, Fm_3, 'LineWidth', 1.7, 'DisplayName', '3 erased');
% plot(par, Fm_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');

legend('Interpreter', 'latex', 'Location','northwest');

xlabel('$T$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$ ','Interpreter', 'latex');

xlim([.5 1])
ylim([.49 1])

% 
% Fm = data.Fm;
% Fmean = mean(transpose(Fm));
% Ncl = sum(transpose(Fm) > F_class)/N;
% 
% % data = load('F_pd_std.mat');
% % Fstd = data.Fstd;
% 
% data = load('F_dir_sigma10.mat');
% Fm_dir = data.Fm_dir;
% Fmean_dir = mean(transpose(Fm_dir));
% Ncl_dir = sum(transpose(Fm_dir) > 0.5)/N;
% 
% % data = load('F_pd_std_dir.mat');
% % Fstd_dir = data.Fstd_dir;
% 
% class = ones(size(par))*F_class;
% 
% set(gca,'fontname','times') 
% subplot(2, 2, [1 3])
% % figure;
% hold on;
% 
% yyaxis left
% % title('Plots with Different y-Scales')
% xlabel('$T$', 'Interpreter', 'latex');
% 
% 
% plot(par, Fmean, 'LineWidth', 1.7, 'DisplayName', 'Protocol');
% plot(par, Fmean_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');
% ylabel('$\bar{\mathcal{F}}$ 'Inter',preter', 'latex');
% 
% yyaxis right
% ylabel('$\%> F_\mathrm{class}$', 'Interpreter', 'latex');
% 
% plot(par, Ncl, 'LineWidth', 1.7, 'DisplayName', 'Protocol');
% plot(par, Ncl_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');
% legend('fontname','times', 'Location','southeast');
% 
% subplot(2, 2, 2)
% title('Protocol')  
% % figure;
% hold on;
% scatter(X, Fm(:, 1:Nsubset), 1, 'MarkerFaceAlpha',.01,'MarkerEdgeAlpha',.01);
% plot(par, Fmean, 'LineWidth', 1.7);
% % plot(x, Ncl, 'LineWidth', 1.7);
% xlabel('$\sigma$', 'Interpreter', 'latex');
% ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
% 
% plot(par, class, 'LineWidth', 1, 'Color', 'red');
% 
% subplot(2, 2,4)
% title('Direct')  
% % figure;
% hold on;
% scatter(X, Fm_dir(:, 1:Nsubset), 1, '+' ,'MarkerFaceAlpha',.01,'MarkerEdgeAlpha',.01);
% plot(par, Fmean_dir, 'LineWidth', 1.7);
% plot(par, class, 'LineWidth', 1, 'Color', 'red');
% 
% % plot(x, Ncl_dir, 'LineWidth', 1.7);
% xlabel('$\sigma$', 'Interpreter', 'latex');
% ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
% 
