clear all;
close all;


fiber_py = load('py_fiber_sig10');
sat_py = load('py_sat_sig10');

tmsv_fiber = fiber_py.data;
tmsv_sat = sat_py.data;

data_fiber = load('results_fiber_sig10');
data_sat = load('results_sat_sig10');

as_fiber = data_fiber.results;
as_sat = data_sat.results;


% Fiber
figure;
hold all;

Ls = linspace(50, 150, 20);
plot(Ls, tmsv_fiber, 'ko-', 'DisplayName', 'TMSV');
plot(Ls, as_fiber, 'o-', 'DisplayName', 'PA-PS');


class = 0.5*ones(length(Ls));
plot(Ls, class, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

legend();
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$L$', 'Interpreter', 'latex');

savefigures('fiber');

% Satellite
figure;
hold all;

zs = [0 5 10 15 20 25 30 35 40 45 50 55 60];

plot(zs, tmsv_sat, 'ko-', 'DisplayName', 'TMSV');
plot(zs, as_sat, 'o-', 'DisplayName', 'PA-PS');

class = 0.5*ones(length(zs));
plot(zs, class, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

legend();
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$\zeta$', 'Interpreter', 'latex');

savefigures('sat');

