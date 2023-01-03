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

par = sqrt(par);


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
ylim([.50 1])

% text(0.96,0.97,'V=3')


savefigures('nonG_V3')
