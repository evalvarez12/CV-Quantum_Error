clear all;
close all;



data = load('results_F_dic.mat');
data = data.results;

Ps = data(:, 1);
F = data(:, 2);
Fdir = data(:, 3);

figure;
hold all;
set(0,'defaultTextInterpreter','latex');


plot(Ps, F, '-', 'LineWidth', 1.7, 'DisplayName', 'EC');
plot(Ps, Fdir, '-','LineWidth', 1.7, 'DisplayName', 'Direct');

legend();
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$P_e$', 'Interpreter', 'latex');

% txt1 = ('$\sigma = 5$');
% text(0.05, .83, txt1);

% ylim([0.45 0.875])



data = load('dic_Fs');
data = data.results;

Vs = data(:, 1);
F1 = data(:, 2);
F2 = data(:, 3);
F3 = data(:, 4);
F12 = data(:, 5);
F13 = data(:, 6);
F23 = data(:, 7);
F123 = F23*0+0.0909;


figure;
hold all;
set(0,'defaultTextInterpreter','latex');

plot(Vs, F1,'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_1$');
plot(Vs, F2,'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_2$');
plot(Vs, F3,'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_3$');
plot(Vs, F12,'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_{12}$');
plot(Vs, F13,'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_{13}$');
plot(Vs, F23,'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_{23}$');
plot(Vs, F23,'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_{23}$');
plot(Vs, F123,'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_{123}$');

l1 = legend();
set(l1, 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$V$', 'Interpreter', 'latex');

