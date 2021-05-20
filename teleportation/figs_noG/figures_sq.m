clear all;
close all;


data_py = load('fixed_py_sq575');
data_py_coh = load('fixed_py');
data_m = load('sq_results_r575');


data_py = data_py.data;
data_m = data_m.results;
data_py_coh = data_py_coh.data;

t = cell2mat(data_py_coh(1));

qs_py = cell2mat(data_py_coh(3));
qs_py_sig2 = qs_py(1,:);
qs_py_sig5 = qs_py(2,:);
qs_py_sig10 = qs_py(3,:);
qs_py_sig20 = qs_py(4,:);

sb_py = data_py;


tmsv_m = data_m(:,1);
ps_m = data_m(:,2);
pa_m = data_m(:,3);
pc_m = data_m(:,4);
as_m = data_m(:,5);
sa_m = data_m(:,6);




figure;
hold all;

% sa_m_sig5(1) = 0.2875;
% ps_m_sig5(1) = 0.27456;
% ps_m_sig5(2) = 0.28067;

plot(t, tmsv_m,'k-', 'LineWidth', 3, 'DisplayName', 'TMSV');
plot(t, sb_py, 'o-', 'DisplayName', 'SB');
plot(t, qs_py_sig5, 'o-', 'DisplayName', 'QS');
plot(t, ps_m, 'o-', 'DisplayName', 'PS');
plot(t, pa_m, 'o-', 'DisplayName', 'PA');
plot(t, pc_m, 'o-', 'DisplayName', 'PC');
plot(t, as_m, 'o-', 'DisplayName', 'PA-PS');
plot(t, sa_m, 'o-', 'DisplayName', 'PS-PA');

class = 0.5*ones(length(t));
plot(t, class, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

legend('Location', 'southeast');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$T$', 'Interpreter', 'latex');

txt1 = ('\bar{s} = 5\mathrm{dB}');
text(0.015, .85, txt1);

% savefigures('fixed1')


