clear all;
close all;



data = load('results_F_dic.mat');
data = data.results;

data_sb = load('results_F_dic_sb.mat');
data_sb = data_sb.results;

Ps = data(:, 1);
F = data(:, 3);
Fdir = data(:, 4);

F_sb = data_sb(:, 3);

figure;
hold all;
set(0,'defaultTextInterpreter','latex');


plot(Ps, F, '-', 'LineWidth', 1.7, 'DisplayName', 'EC');
plot(Ps, F, 'o', 'LineWidth', 1.7, 'DisplayName', 'SB');
plot(Ps, Fdir, '-','LineWidth', 1.7, 'DisplayName', 'Direct');

legend();
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$P_e$', 'Interpreter', 'latex');

% txt1 = ('$\sigma = 5$');
% text(0.05, .83, txt1);

% ylim([0.45 0.875])



data = load('dic_Fs');
data = data.results;

data_sb = load('dic_Fs_sb');
data_sb = data_sb.results;

Vs = data(:, 1);
F1 = data(:, 2);
F2 = data(:, 3);
F3 = data(:, 4);
F12 = data(:, 5);
F13 = data(:, 6);
F23 = data(:, 7);
F123 = F23*0+0.0909;


F1_sb = data_sb(:, 2);
F2_sb = data_sb(:, 3);
F3_sb = data_sb(:, 4);
F12_sb = data_sb(:, 5);
F13_sb = data_sb(:, 6);
F23_sb = data_sb(:, 7);


figure;
hold all;
set(0,'defaultTextInterpreter','latex');

plot(Vs(1:2:end), F1(1:2:end),'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_1$');
plot(Vs(1:2:end), F2(1:2:end),'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_2$');
plot(Vs(1:2:end), F3(1:2:end),'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_3$');
plot(Vs(1:2:end), F12(1:2:end),'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_{12}$');
plot(Vs(1:2:end), F13(1:2:end),'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_{13}$');
plot(Vs(1:2:end), F23(1:2:end),'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_{23}$');
plot(Vs(1:2:end), F23(1:2:end),'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_{23}$');
plot(Vs(1:2:end), F123(1:2:end),'-', 'LineWidth', 1.7, 'DisplayName', '$\bar{\mathcal{F}}_{123}$');


plot(Vs(1:2:end), F1_sb(1:2:end),'o', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB}~ \bar{\mathcal{F}}_1$');
plot(Vs(1:2:end), F2_sb(1:2:end),'o', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB}~ \bar{\mathcal{F}}_2$');
plot(Vs(1:2:end), F3_sb(1:2:end),'o', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB}~\bar{\mathcal{F}}_3$');
plot(Vs(1:2:end), F12_sb(1:2:end),'o', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB}~\bar{\mathcal{F}}_{12}$');
plot(Vs(1:2:end), F13_sb(1:2:end),'o', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB}~\bar{\mathcal{F}}_{13}$');
plot(Vs(1:2:end), F23_sb(1:2:end),'o', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB}~\bar{\mathcal{F}}_{23}$');
plot(Vs(1:2:end), F23_sb(1:2:end),'o', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB}~\bar{\mathcal{F}}_{23}$');
plot(Vs(1:2:end), F123(1:2:end),'o', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB}~\bar{\mathcal{F}}_{123}$');

l1 = legend('NumColumns',4);
set(l1, 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$V$', 'Interpreter', 'latex');

