clear all;
close all;
set(0,'defaultTextInterpreter','latex');

%par = linspace(0.4, 1, 24);
% Data structure
% results = [Ts(:), F_tmsv(:), F_ps(:), F_pa(:), F_pc(:), F_pspa(:), F_paps(:)];


sigma_coh = 10;
F_class = (1 + (1/sigma_coh))/(2 + (1/sigma_coh));

% dataV32 = load('fid_nonG_T1erased_V3_s10.mat');

dataV3 = load('fid_T1erased_V3_s10.mat');
dataV6 = load('fid_T1erased_V6_s10.mat');
dataV10 = load('fid_T1erased_V10_s10.mat');
% dataV20 = load('fid_nonG_T2erased_V20_s10.mat');

dataV3 = dataV3.results;
dataV6 = dataV6.results;
dataV10 = dataV10.results;
% dataV20 = dataV20.results;

data = dataV3;

par = data(:, 1);
tmsv = data(:, 2);
ps = data(:, 3);
pa = data(:, 4);
pc = data(:, 5);
pspa = data(:, 6);
paps = data(:, 7);
sb = data(:, 8);


hold on;
plot(par, ps, 'o-', 'LineWidth', 1.7, 'DisplayName', 'PS');
plot(par, pa, '+-','LineWidth', 1.7, 'DisplayName', 'PA');
plot(par, pc, 'v-', 'LineWidth', 1.7, 'DisplayName', 'PC');
plot(par, pspa, '*-', 'LineWidth', 1.7, 'DisplayName', 'PS-PA');
plot(par, paps, '.-', 'LineWidth', 1.7, 'DisplayName', 'PA-PS');
plot(par, sb, '^-',  'LineWidth', 1.7, 'DisplayName', 'SB');
plot(par, tmsv, 'k-',  'LineWidth', 3, 'DisplayName', 'TMSV');

legend('Location','northwest');

xlabel('$T$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$ ','Interpreter', 'latex');

% xlim([.5 1])
ylim([.30 1])

text(0.96,0.97,'V=3')


savefigures('nonG_V3')

% 
% F1V1 = dataV1(:, :, 2);
% Fm_1V1= mean(F1V1, 2);
% 
% 
% dataV3 = load('F_compare_V3_sigma10_sigma10.mat');
% % dataV3 = load('F_compare_V3_sigma10_sigma0.mat');
% % dataV3 = load('F_compare_V3_sigma10_sigma05.mat');
% 
% dataV3 = dataV3.results;
% 
% FnoV3 = dataV3(:, :, 1);
% F1V3 = dataV3(:, :, 2);
% F2V3 = dataV3(:, :, 3);
% F3V3 = dataV3(:, :, 4);
% FdirV3= dataV3(:, :, 5);
% 
% Fm_noV3 = mean(FnoV3, 2);
% Fm_1V3= mean(F1V3, 2);
% Fm_2V3= mean(F2V3, 2);
% Fm_3V3= mean(F3V3, 2);
% Fm_dirV3 =  mean(FdirV3, 2);
% 
% 
% dataV5 = load('F_compare_V5_sigma10_sigma10.mat');
% % dataV5 = load('F_compare_V5_sigma10_sigma0.mat');
% % dataV5 = load('F_compare_V5_sigma10_sigma05.mat');
% 
% dataV5 = dataV5.results;
% 
% F1V5 = dataV5(:, :, 2);
% Fm_1V5= mean(F1V5, 2);
% 
% dataV10 = load('F_compare_V10_sigma10_sigma10.mat');
% % dataV10 = load('F_compare_V10_sigma10_sigma0.mat');
% % dataV10 = load('F_compare_V10_sigma10_sigma05.mat');
% 
% dataV10 = dataV10.results;
% 
% F1V10 = dataV10(:, :, 2);
% Fm_1V10= mean(F1V10, 2);
% 
% dataV20 = load('F_compare_V20_sigma10_sigma10.mat');
% % dataV20 = load('F_compare_V20_sigma10_sigma0.mat');
% % dataV20 = load('F_compare_V20_sigma10_sigma05.mat');
% 
% dataV20 = dataV20.results;
% 
% F1V20 = dataV20(:, :, 2);
% Fm_1V20= mean(F1V20, 2);
% 
% FnoV20 = dataV20(:, :, 1);
% Fm_noV20 = mean(FnoV20, 2);
% 
% Fdir = dataV20(:, :, 5);
% Fm_dir = mean(Fdir, 2);
% 
% 
% set(gca,'fontname','times') 
% hold on;
% 
% plot(par, Fm_noV3, 'LineWidth', 1.7, 'DisplayName', 'No erasure');
% plot(par, Fm_1V3, 'LineWidth', 1.7, 'DisplayName', 'Erasure - $V=3$');
% plot(par, Fm_1V5, 'LineWidth', 1.7, 'DisplayName', 'Erasure - $V=5$');
% plot(par, Fm_1V10, 'LineWidth', 1.7, 'DisplayName', 'Erasure - $V=10$');
% plot(par, Fm_1V20, 'LineWidth', 1.7, 'DisplayName', 'Erasure - $V=20$');
% 
% 
% legend('Interpreter', 'latex', 'Location','northwest');
% 
% xlabel('$T$', 'Interpreter', 'latex');
% ylabel('$\bar{\mathcal{F}}$ ','Interpreter', 'latex');
% 
% xlim([.5 1])
% ylim([.49 1])
