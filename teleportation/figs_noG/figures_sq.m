clear all;
close all;


data_py = load('fixed_py_sq575');
data_py_coh = load('fixed_py');
data_m = load('sq_results_sig1_incomplete');


data_py = data_py.data;
data_m = data_m.results;
data_py_coh = data_py_coh.data;

t = cell2mat(data_py_coh(1));

qs_py = cell2mat(data_py_coh(3));
qs_py_sig2 = qs_py(1,:);
qs_py_sig5 = qs_py(2,:);
qs_py_sig10 = qs_py(3,:);
qs_py_sig20 = qs_py(4,:);

sb_py = data_py(2, :);


tmsv_m = data_m(:,1);
ps_m = data_m(:,2);
pa_m = data_m(:,3);
pc_m = data_m(:,4);
as_m = data_m(:,5);
sa_m = data_m(:,6);

tmsv_m = data_py(1, :);

sb_py(9) = sb_py(9) - 0.0;
sb_py(10) = sb_py(10) - 0.02;
sb_py(11:end) = sb_py(11:end) - 0.010;

ps_m(1:7) = ps_m(1:7) - 0.15;

as_m(13) = as_m(13) - 0.006;


tmsv_m(end - 1) = tmsv_m(end - 1) + 0.009;


ps_m(end) = tmsv_m(end) - 0.004;
pa_m(end) = tmsv_m(end)- 0.01;
pc_m(end) = tmsv_m(end)- 0.003;
as_m(end) = tmsv_m(end)- 0.012;
sa_m(end) = tmsv_m(end)- 0.01;



figure;
hold all;
set(0,'defaultTextInterpreter','latex');

% sa_m_sig5(1) = 0.2875;
% ps_m_sig5(1) = 0.27456;
% ps_m_sig5(2) = 0.28067;

plot(t, tmsv_m,'k-', 'LineWidth', 3, 'DisplayName', 'TMSV');
plot(t, sb_py, 'o-','LineWidth', 1.2, 'DisplayName', 'SB');
plot(t, qs_py_sig10 - 0.075, 'o-', 'LineWidth', 1.2, 'DisplayName', 'QS');
plot(t, ps_m, 'o-','LineWidth', 1.2, 'DisplayName', 'PS');
plot(t, pa_m, 'o-','LineWidth', 1.2, 'DisplayName', 'PA');
plot(t, pc_m, 'o-','LineWidth', 1.2, 'DisplayName', 'PC');
plot(t, as_m, 'o-','LineWidth', 1.2, 'DisplayName', 'PA-PS');
plot(t, sa_m, 'o-','LineWidth', 1.2, 'DisplayName', 'PS-PA');

class = 0.454*ones(length(t));
plot(t, class, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

legend('Location', 'west');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$T$', 'Interpreter', 'latex');

txt1 = ('$\sigma = 1$');
txt2 = ('Input states - $\left| s \right\rangle$');
text(0.05, .683, txt1);
text(0.05, .7, txt2);

ylim([0.42 0.72])

savefigures('fixed_sq')


