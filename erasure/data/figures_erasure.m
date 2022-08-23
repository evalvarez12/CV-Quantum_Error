clear all;
close all;
set(0,'defaultTextInterpreter','latex');

data1 = load('F_erasure_T1.mat');
data1 = data1.Fm;

data2 = load('F_erasure_T2.mat');
data2 = data2.Fm;

data3 = load('F_erasure_T3.mat');
data3 = data3.Fm;

datadir = load('F_erasure_dir.mat');
datadir = datadir.Fm_dir;

x = 0.3:0.05:1;

sigma_coh = 10;
F_class = (1 + (1/sigma_coh))/(2 + (1/sigma_coh));



set(gca,'fontname','times') 
hold on;
plot(x, data1, 'LineWidth', 1.2, 'DisplayName', 'T1');
plot(x, data2,'o-', 'LineWidth', 1.7, 'DisplayName', 'T2');
plot(x, data3,'+-', 'LineWidth', 1.7, 'DisplayName', 'T3');

plot(x, datadir, 'LineWidth', 1.7, 'DisplayName', 'dir');

xlabel('$T$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');

