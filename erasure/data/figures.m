clear all;
close all;


% define colors
red1= '#D10303';
red2= '#F45F50';

blue1 = '#071AB3';
blue2 = '#4F6FFC';


data = load('results_F_dic.mat');
data = data.results;

data_sb = load('results_F_dic_sb.mat');
data_sb = data_sb.results;

Ps = data(:, 1);
F = data(:, 3);
Fdir = data(:, 4);

F_sb = data_sb(:, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fidelities plot

figure;
hold all;
set(0,'defaultTextInterpreter','latex');


plot(Ps, F, '--', 'LineWidth', 1.7, 'DisplayName', 'TMSV', 'Color', blue1);
plot(Ps, F_sb, 'o-', 'LineWidth', 1.7, 'DisplayName', 'SB', 'Color', red1);
plot(Ps, Fdir, '-','LineWidth', 1.7, 'DisplayName', 'Direct', 'Color', 'black');

l1 = legend();
set(l1, 'Interpreter', 'latex');
ylabel('$\mathcal{F}_\mathrm{total}$', 'Interpreter', 'latex');
xlabel('$P_e$', 'Interpreter', 'latex');

% txt1 = ('$\sigma = 5$');
% text(0.05, .83, txt1);

xlim([0. 0.5])



%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal squeezing plot

V = data(:, 2);
V_sb = data_sb(:, 2);

V_sb(4) = V_sb(3);

% figure;
% hold all;
% set(0,'defaultTextInterpreter','latex');

axes('Position',[.17 .18 .4 .4])
box on

hold all;
plot(Ps, V, '--', 'LineWidth', 1.7, 'DisplayName', 'TMSV', 'Color', blue1);
plot(Ps, V_sb, 'o-', 'LineWidth', 1.7, 'DisplayName', 'SB', 'Color', red1);

legend();
ylabel('$V_\mathrm{opt}$', 'Interpreter', 'latex');
xlabel('$P_e$', 'Interpreter', 'latex');
xlim([0. 0.5])


% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fidelity components plot
% 
% data = load('dic_Fs');
% data = data.results;
% 
% data_sb = load('dic_Fs_sb');
% data_sb = data_sb.results;
% 
% Vs = data(:, 1);
% F1 = data(:, 2);
% F2 = data(:, 3);
% F3 = data(:, 4);
% F12 = data(:, 5);
% F13 = data(:, 6);
% F23 = data(:, 7);
% F123 = F23*0+0.0909;
% 
% 
% F1_sb = data_sb(:, 2);
% F2_sb = data_sb(:, 3);
% F3_sb = data_sb(:, 4);
% F12_sb = data_sb(:, 5);
% F13_sb = data_sb(:, 6);
% F23_sb = data_sb(:, 7);
% 
% 
% % Fixing deviations
% F1_sb(50:72) = F2_sb(50:72);
% 
% figure;
% hold all;
% set(0,'defaultTextInterpreter','latex');
% 
% plot(Vs(1:3:end), F1(1:3:end),'-', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{TMSV-}\mathcal{F}_1$', 'Color', blue1);
% plot(Vs(1:3:end), F2(1:3:end),'o', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{TMSV-}\mathcal{F}_2$', 'Color', blue2);
% % plot(Vs(1:2:end), F3(1:2:end),'-', 'LineWidth', 1.7, 'DisplayName', '$3$');
% % plot(Vs(1:2:end), F12(1:2:end),'o', 'LineWidth', 1.7, 'DisplayName', '${12}$');
% plot(Vs(1:3:end), F13(1:3:end),'d', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{TMSV-}\mathcal{F}_{13}$', 'Color', blue1);
% plot(Vs(1:3:end), F23(1:3:end),'p', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{TMSV-}\mathcal{F}_{23}$', 'Color', blue2);
% % plot(Vs(1:2:end), F123(1:2:end),'-', 'LineWidth', 1.7, 'DisplayName', '${123}$');
% 
% 
% plot(Vs(1:3:end), F1_sb(1:3:end),'--', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}\mathcal{F}_1$', 'Color', red1);
% plot(Vs(1:3:end), F2_sb(1:3:end),'s', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}\mathcal{F}_2$', 'Color', red2);
% % plot(Vs(1:2:end), F3_sb(1:2:end),'--', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}3$');
% % plot(Vs(1:2:end), F12_sb(1:2:end),'+', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}{12}$');
% plot(Vs(1:3:end), F13_sb(1:3:end),'v', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}\mathcal{F}_{13}$', 'Color', red1);
% plot(Vs(1:3:end), F23_sb(1:3:end),'x', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}\mathcal{F}_{23}$', 'Color', red2);
% % plot(Vs(1:2:end), F123(1:2:end),'--', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}{123}$');
% 
% l1 = legend('NumColumns',2);
% set(l1, 'Interpreter', 'latex');
% ylabel('$\mathcal{F}$', 'Interpreter', 'latex');
% xlabel('$V$', 'Interpreter', 'latex');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Optimal g component plot
% 
% data = load('dic_Fs_pars2');
% data = data.results_pars2;
% 
% data_sb = load('dic_Fs_pars_sb');
% data_sb = data_sb.results_pars;
% 
% 
% F1 = data(:, 1);
% F2 = data(:, 2);
% F3 = data(:, 3);
% F12 = data(:, 4);
% F13 = data(:, 5);
% F23 = data(:, 6);
% 
% F1_sb = data_sb(:, 1);
% F2_sb = data_sb(:, 3);
% F3_sb = data_sb(:, 5);
% F12_sb = data_sb(:, 7);
% F13_sb = data_sb(:, 9);
% F23_sb = data_sb(:, 11);
% 
% % Fixing deviations
% F1_sb(59) = -F2_sb(59);
% F1_sb(61) = -F2_sb(61);
% F1_sb(67) = -F2_sb(67);
% 
% 
% figure;
% hold all;
% set(0,'defaultTextInterpreter','latex');
% 
% plot(Vs(1:2:end), F1(1:2:end), '-','LineWidth', 1.7, 'DisplayName', '$\mathrm{TMSV-}\mathcal{F}_1$', 'Color', blue1);
% plot(Vs(1:3:end), F2(1:3:end), 'o', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{TMSV-}\mathcal{F}_2$', 'Color', blue2);
% % plot(Vs(1:2:end), F3(1:2:end), 'LineWidth', 1.7, 'DisplayName', '$3$');
% % plot(Vs(1:2:end), F12(1:2:end), 'LineWidth', 1.7, 'DisplayName', '${12}$');
% plot(Vs(1:2:end), F13(1:2:end), 'd','LineWidth', 1.7, 'DisplayName', '$\mathrm{TMSV-}\mathcal{F}_{13}$', 'Color', blue1);
% plot(Vs(1:3:end), F23(1:3:end), 'p', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{TMSV-}\mathcal{F}_{23}$', 'Color', blue2);
% % plot(Vs(1:2:end), F123(1:2:end), 'LineWidth', 1.7, 'DisplayName', '${123}$');
% 
% 
% plot(Vs(1:2:end), F1_sb(1:2:end), '--', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}\mathcal{F}_1$', 'Color', red1);
% plot(Vs(1:3:end), F2_sb(1:3:end), 's', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}\mathcal{F}_2$', 'Color', red2);
% % plot(Vs(1:2:end), F3_sb(1:2:end), '--', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}3$');
% % plot(Vs(1:2:end), F12_sb(1:2:end), '+', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}{12}$');
% plot(Vs(1:2:end), F13_sb(1:2:end), 'v', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}\mathcal{F}_{13}$', 'Color', red1);
% plot(Vs(1:3:end), F23_sb(1:3:end), 'x', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}\mathcal{F}_{23}$', 'Color', red2);
% % plot(Vs(1:2:end), F123(1:2:end), '--', 'LineWidth', 1.7, 'DisplayName', '$\mathrm{SB-}{123}$');
% 
% l1 = legend('NumColumns',2);
% set(l1, 'Interpreter', 'latex');
% ylabel('$g_\mathrm{opt}$', 'Interpreter', 'latex');
% xlabel('$V$', 'Interpreter', 'latex');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Postselected fidelities plot
% 
% 
% %%%Component values
% %F1 V=2 0.8295
% %F1 V=4 0.9155
% %F1 V=6 0.9459
% 
% %SB F1 V=2 0.8988
% %SB F1 V=4 0.9532
% %SB F1 V=6 0.9711
% 
% data9 = load('results_F_dic_part.mat');
% data6 = load('results_F_dic_part_v6.mat');
% data3 = load('results_F_dic_part_v3.mat');
% data9 = data9.results;
% data6 = data6.results;
% data3 = data3.results;
% 
% 
% data9sb = load('results_F_dic_part_v9_sb.mat');
% data6sb = load('results_F_dic_part_v6_sb.mat');
% data3sb = load('results_F_dic_part_v3_sb.mat');
% data9sb = data9sb.results;
% data6sb = data6sb.results;
% data3sb = data3sb.results;
% 
% Ps = data9(:, 1);
% Fdir = data9(:, 4);
% 
% Pnorm = (1-Ps).^3 + 3.*Ps.*(1-Ps).^2; 
% % Pnorm = Ps.^3 + 3.*Ps.^2.*(1-Ps); 
% 
% % F9 = data9(:, 3) ./ Pnorm;
% % F6 = -data6(:,3) ./ Pnorm;
% % F3 = -data3(:, 3) ./ Pnorm;
% 
% 
% 
% F9 = data9(:, 3);
% F6 = -data6(:,3);
% F3 = -data3(:, 3);
% 
% F9sb = -data9sb(:, 3);
% F6sb = -data6sb(:,3);
% F3sb = -data3sb(:, 3);
% 
% 
% % F9 = (1-Ps) + Ps.*(2/3 * 0.9459 + 1/3);
% % F6 = (1-Ps) + Ps.*(2/3 * 0.9155 + 1/3);
% % F3 = (1-Ps) + Ps.*(2/3 * 0.8295 + 1/3);
% % 
% % F9sb = (1-Ps) + Ps.*(2/3 * 0.9711 + 1/3);
% % F6sb = (1-Ps) + Ps.*(2/3 * 0.9532 + 1/3);
% % F3sb = (1-Ps) + Ps.*(2/3 * 0.8988 + 1/3);
% 
% figure;
% hold all;
% set(0,'defaultTextInterpreter','latex');
% 
% 
% plot(Ps, F, '--', 'LineWidth', 1.7, 'DisplayName', 'Optimized');
% 
% plot(Ps, F9, 'o-', 'LineWidth', 1.7, 'DisplayName', '$V=9$');
% plot(Ps, F6, 's-', 'LineWidth', 1.7, 'DisplayName', '$V=6$');
% plot(Ps, F3, 'd-', 'LineWidth', 1.7, 'DisplayName', '$V=3$');
% 
% % 
% % plot(Ps, F9sb, 'o-', 'LineWidth', 1.7, 'DisplayName', 'SB');
% % plot(Ps, F6sb, '+-', 'LineWidth', 1.7, 'DisplayName', 'SB');
% % plot(Ps, F3sb, 'v-', 'LineWidth', 1.7, 'DisplayName', 'SB');
% 
% 
% plot(Ps, Fdir, '-','LineWidth', 1.7, 'DisplayName', 'Direct', 'Color', 'black');
% 
% l1 = legend();
% set(l1, 'Interpreter', 'latex');
% ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
% xlabel('$P_e$', 'Interpreter', 'latex');
% 
% % txt1 = ('$\sigma = 5$');
% % text(0.05, .83, txt1);
% 
% % ylim([0.45 0.875])
