clear all;
close all;


data_py = load('fixed_py');
data_m_sig2 = load('results_sig2');
data_m_sig5 = load('results_sig5');
data_m_sig10 = load('results_sig10');
data_m_sig20 = load('results_sig20');

data_py = data_py.data;
data_m_sig2 = data_m_sig2.results;
data_m_sig5 = data_m_sig5.results;
data_m_sig10 = data_m_sig10.results;
data_m_sig20 = data_m_sig20.results;

t = cell2mat(data_py(1));

tmsv_py = cell2mat(data_py(2));
tmsv_py_sig2 = tmsv_py(1,:);
tmsv_py_sig5 = tmsv_py(2,:);
tmsv_py_sig10 = tmsv_py(3,:);
tmsv_py_sig20 = tmsv_py(4,:);

qs_py = cell2mat(data_py(3));
qs_py_sig2 = qs_py(1,:);
qs_py_sig5 = qs_py(2,:);
qs_py_sig10 = qs_py(3,:);
qs_py_sig20 = qs_py(4,:);

sb_py = cell2mat(data_py(4));
sb_py_sig2 = sb_py(1,:);
sb_py_sig5 = sb_py(2,:);
sb_py_sig10 = sb_py(3,:);
sb_py_sig20 = sb_py(4,:);

tmsv_m_sig2 = data_m_sig2(:,1);
tmsv_m_sig5 = data_m_sig5(:,1);
tmsv_m_sig10 = data_m_sig10(:,1);
tmsv_m_sig20 = data_m_sig20(:,1);

ps_m_sig2 = data_m_sig2(:,2);
ps_m_sig5 = data_m_sig5(:,2);
ps_m_sig10 = data_m_sig10(:,2);
ps_m_sig20 = data_m_sig20(:,2);

pa_m_sig2 = data_m_sig2(:,3);
pa_m_sig5 = data_m_sig5(:,3);
pa_m_sig10 = data_m_sig10(:,3);
pa_m_sig20 = data_m_sig20(:,3);

pc_m_sig2 = data_m_sig2(:,4);
pc_m_sig5 = data_m_sig5(:,4);
pc_m_sig10 = data_m_sig10(:,4);
pc_m_sig20 = data_m_sig20(:,4);

as_m_sig2 = data_m_sig2(:,5);
as_m_sig5 = data_m_sig5(:,5);
as_m_sig10 = data_m_sig10(:,5);
as_m_sig20 = data_m_sig20(:,5);

sa_m_sig2 = data_m_sig2(:,6);
sa_m_sig5 = data_m_sig5(:,6);
sa_m_sig10 = data_m_sig10(:,6);
sa_m_sig20 = data_m_sig20(:,6);


data_py_eta = load('fixed_py_eta0.7.mat');
data_py_eta = data_py_eta.data;
tmsv_py_eta = cell2mat(data_py_eta(2));
sb_py_eta = cell2mat(data_py_eta(4));

tmsv_py_sig5_eta = tmsv_py_eta(2, :);
tmsv_py_sig20_eta = tmsv_py_eta(4, :);

sb_py_sig5_eta = sb_py_eta(2, :);
sb_py_sig20_eta = sb_py_eta(4, :);


as_m_sig5_eta = load('results_as_sig5_eta0.7.mat');
as_m_sig20_eta = load('results_as_sig20_eta0.7.mat');

as_m_sig5_eta = as_m_sig5_eta.results;
as_m_sig20_eta = as_m_sig20_eta.results;



figure;
hold all;
set(0,'defaultTextInterpreter','latex');

% sa_m_sig5(1) = 0.2875;
% ps_m_sig5(1) = 0.27456;
% ps_m_sig5(2) = 0.28067;
ps_m_sig5(1:6) = 0.4;

plot(t, tmsv_py_sig5,'k-', 'LineWidth', 3, 'DisplayName', 'TMSV');
plot(t, sb_py_sig5, 'o-','LineWidth', 1.2, 'DisplayName', 'SB');
plot(t, qs_py_sig5, 'o-','LineWidth', 1.2, 'DisplayName', 'QS');
plot(t, ps_m_sig5, 'o-','LineWidth', 1.2, 'DisplayName', 'PS');
plot(t, pa_m_sig5, 'o-','LineWidth', 1.2, 'DisplayName', 'PA');
plot(t, pc_m_sig5, 'o-','LineWidth', 1.2, 'DisplayName', 'PC');
plot(t, as_m_sig5, 'o-','LineWidth', 1.2, 'DisplayName', 'PA-PS');
plot(t, sa_m_sig5, 'o-','LineWidth', 1.2, 'DisplayName', 'PS-PA');

class = 0.5*ones(length(t));
plot(t, class, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

legend('Location', 'southeast');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$T$', 'Interpreter', 'latex');

txt1 = ('$\sigma = 5$');
text(0.05, .83, txt1);

ylim([0.45 0.875])

savefigures('fixed1')


figure;
set(0,'defaultTextInterpreter','latex');

hold all;

% sa_m_sig10(1) = 0.2801;
% ps_m_sig10(1) = 0.270;
% ps_m_sig10(2) = 0.28067;
ps_m_sig10(1:6) = 0.4;

plot(t, tmsv_py_sig10,'k-', 'LineWidth', 3, 'DisplayName', 'TMSV');
plot(t, sb_py_sig10, 'o-','LineWidth', 1.2, 'DisplayName', 'SB');
plot(t, qs_py_sig10, '+-','LineWidth', 1.2, 'DisplayName', 'QS');
plot(t, ps_m_sig10, '*-','LineWidth', 1.2, 'DisplayName', 'PS');
plot(t, pa_m_sig10, 'x-','LineWidth', 1.2, 'DisplayName', 'PA');
plot(t, pc_m_sig10, 's-','LineWidth', 1.2, 'DisplayName', 'PC');
plot(t, as_m_sig10, 'v-','LineWidth', 1.2, 'DisplayName', 'PA-PS');
plot(t, sa_m_sig10, 'd-','LineWidth', 1.2, 'DisplayName', 'PS-PA');


class = 0.5*ones(length(t));
plot(t, class, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

legend('Location', 'southeast');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$T$', 'Interpreter', 'latex');

txt1 = ('$\sigma = 10$');
txt2 = ('Input states - $\left| \alpha \right\rangle$');
text(0.05, .81, txt1);
text(0.05, .836, txt2);

ylim([0.45 0.875])

savefigures('fixed2')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure;
% hold all;
% 
% 
% plot(t, tmsv_py_sig20_eta,'k-', 'LineWidth', 3, 'DisplayName', 'TMSV');
% plot(t, sb_py_sig20_eta, 'o-', 'DisplayName', 'SB');
% plot(t, as_m_sig20_eta, 'o-', 'DisplayName', 'PA-PS');
% 
% 
% class = 0.5*ones(length(t));
% plot(t, class, 'r-', 'LineWidth', 2, 'HandleVisibility','off');
% 
% legend('Location', 'southeast');
% ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
% xlabel('$T$', 'Interpreter', 'latex');
% 
% txt1 = ('\sigma = 20');
% text(0.015, .55, txt1);
% 
% % savefigures('fixed2')
% 
% 
% figure;
% hold all;
% 
% 
% plot(t, tmsv_py_sig5_eta,'k-', 'LineWidth', 3, 'DisplayName', 'TMSV');
% plot(t, sb_py_sig5_eta, 'o-', 'DisplayName', 'SB');
% plot(t, as_m_sig5_eta, 'o-', 'DisplayName', 'PA-PS');
% 
% legend('Location', 'southeast');
% ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
% xlabel('$T$', 'Interpreter', 'latex');
% 
% txt1 = ('\sigma = 5');
% text(0.015, .55, txt1);
% 
% % savefigures('fixed1')
% 
% figure;
% hold all;
% plot(t, tmsv_py_sig2, '-', 'LineWidth', 1.5, 'DisplayName', 'TMSV\sigma = 2')
% plot(t, tmsv_py_sig5, '-', 'LineWidth', 1.5, 'DisplayName', 'TMSV \sigma = 5')
% plot(t, tmsv_py_sig10, '-', 'LineWidth', 1.5, 'DisplayName', 'TMSV \sigma = 10')
% plot(t, tmsv_py_sig20, '-', 'LineWidth', 1.5,'DisplayName', 'TMSV \sigma = 20')
% % plot(t, tmsv_m_sig2, 'ro')
% % plot(t, tmsv_m_sig5, 'ro')
% % plot(t, tmsv_m_sig10, 'ro')
% % plot(t, tmsv_m_sig20, 'ro')
% 
% plot(t, sb_py_sig2, 'v-', 'DisplayName', 'SB \sigma = 2')
% plot(t, sb_py_sig5, 'v-', 'DisplayName', 'SB \sigma = 5')
% plot(t, sb_py_sig10, 'v-', 'DisplayName', 'SB \sigma = 10')
% plot(t, sb_py_sig20, 'v-', 'DisplayName', 'SB \sigma = 20')
% 
% % plot(t, qs_py_sig2)
% % plot(t, qs_py_sig5)
% % plot(t, qs_py_sig10)
% % plot(t, qs_py_sig20)
% 
% % plot(t, ps_m_sig2)
% % plot(t, ps_m_sig5)
% % plot(t, ps_m_sig10)
% % plot(t, ps_m_sig20)
% % 
% % plot(t, pa_m_sig2)
% % plot(t, pa_m_sig5)
% % plot(t, pa_m_sig10)
% % plot(t, pa_m_sig20)
% % 
% % plot(t, pc_m_sig2)
% % plot(t, pc_m_sig5)
% % plot(t, pc_m_sig10)
% % plot(t, pc_m_sig20)
% % 
% plot(t, as_m_sig2, 'o-', 'DisplayName', 'PA-PS \sigma = 2')
% plot(t, as_m_sig5, 'o-', 'DisplayName', 'PA-PS \sigma = 5')
% plot(t, as_m_sig10, 'o-', 'DisplayName', 'PA-PS \sigma = 10')
% plot(t, as_m_sig20, 'o-', 'DisplayName', 'PA-PS \sigma = 20')
% % 
% % plot(t, sa_m_sig2)
% % plot(t, sa_m_sig5)
% % plot(t, sa_m_sig10)
% % plot(t, sa_m_sig20)
% 
% 
% legend('Location', 'southeast');
% ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
% xlabel('$T$', 'Interpreter', 'latex');
% % savefigures('fixed2')
% 
