clear all;
close all;
set(0,'defaultTextInterpreter','latex');

vs = 1.5e12;
eps = 0.013;
eta = 1;
T = linspace(0.9999999, .6, 100);
V = 1.5e12;

pe = 0.;
R0 = RCI(V, vs, T, eps, eta, pe);
Rdir0 = RCI_dir(vs, T, eps, pe);

pe = 0.1;
R1 = RCI(V, vs, T, eps, eta, pe);
Rdir1 = RCI_dir(vs, T, eps, pe);

pe = 0.2;
R2 = RCI(V, vs, T, eps, eta, pe);
Rdir2 = RCI_dir(vs, T, eps, pe);

pe = 0.3;
R3 = RCI(V, vs, T, eps, eta, pe);
Rdir3 = RCI_dir(vs, T, eps, pe);

color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];
color3 = [0.9290, 0.6940, 0.1250];

hold on;
plot(T, R0, '-', 'Color', 'black', 'LineWidth', 1.6, 'DisplayName', '$p_\mathrm{e}=0.$');
plot(T, Rdir0,'--', 'Color', 'black', 'LineWidth', 1.4, 'HandleVisibility','off');

plot(T, R1, '-', 'Color', color1, 'LineWidth', 1.4, 'DisplayName', '$p_\mathrm{e}=0.1$');
plot(T, Rdir1,'--', 'Color', color1, 'LineWidth', 1.4, 'HandleVisibility','off');

plot(T, R2, '-', 'Color', color2, 'LineWidth', 1.4, 'DisplayName', '$p_\mathrm{e}=0.2$');
plot(T, Rdir2,'--', 'Color', color2, 'LineWidth', 1.4, 'HandleVisibility','off');

plot(T, R3, '-', 'Color', color3, 'LineWidth', 1.4, 'DisplayName', '$p_\mathrm{e}=0.3$');
plot(T, Rdir3,'--', 'Color', color3, 'LineWidth', 1.4, 'HandleVisibility','off');

set(gca, 'YScale', 'log');
set(gca,'fontname','times') 

xlabel('$T$', 'Interpreter', 'latex');
ylabel('$\mathcal{R}$', 'Interpreter', 'latex');
legend('fontname','times', 'Location','northwest','Interpreter', 'latex');

xlim([0.8 1]);
% ylim([1 10]);
savefigures('RCI')

% 
% 
% set(gca,'fontname','times') 
% 
% figure;
% hold on;
% plot(par_erasure, data3_erasure, 'p-', 'LineWidth', 1.5, 'DisplayName', 'V=10');
% plot(par_erasure, data2_erasure, '*-', 'LineWidth', 1.5, 'DisplayName', 'V=6');
% plot(par_erasure, data1_erasure, 'v-', 'LineWidth', 1.5, 'DisplayName', 'V=3');
% 
% 
% % plot(par_no, data1_no, 'LineWidth', 1.2, 'DisplayName', 'V=3');
% plot(par_no, data2_no, 'LineWidth', 1.5, 'DisplayName', 'No erasure');
% 
% ylim([0.4 1]);
% 
% xlabel('$T$', 'Interpreter', 'latex');
% ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
% % legend('Location','northwest')
% text(0.95, 0.97, '$\sigma=0.2$')
% 
% savefigures('erasure_sigma2')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% data1_erasure = load('Ferasure_T1_V3_sigma0.mat');
% data1_erasure = transpose(data1_erasure.Fm);
% 
% data2_erasure = load('Ferasure_T1_V6_sigma0.mat');
% data2_erasure = transpose(data2_erasure.Fm);
% 
% 
% data3_erasure = load('Ferasure_T1_V10_sigma0.mat');
% data3_erasure = transpose(data3_erasure.Fm);
% 
% data1_no = load('Fnoerasure_V3_sigma0.mat');
% data1_no = transpose(data1_no.Fm);
% 
% 
% 
% 
% figure;
% hold on;
% plot(par_erasure, data3_erasure, 'p-', 'LineWidth', 1.5, 'DisplayName', 'V=10');
% plot(par_erasure, data2_erasure, '*-', 'LineWidth', 1.5, 'DisplayName', 'V=6');
% plot(par_erasure, data1_erasure, 'v-', 'LineWidth', 1.5, 'DisplayName', 'V=3');
% 
% 
% % plot(par_no, data1_no, 'LineWidth', 1.2, 'DisplayName', 'V=3');
% plot(par_erasure, data1_no, 'LineWidth', 1.5, 'DisplayName', 'No erasure');
% 
% 
% 
% xlabel('$T$', 'Interpreter', 'latex');
% ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
% legend('Location','northwest')
% text(0.81, 0.97, '$\sigma=0$')
% 
% 
% savefigures('erasure_sigma0')
