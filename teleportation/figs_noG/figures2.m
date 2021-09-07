clear all;
close all;


fiber_py = load('py_fiber_sig10-2');
sat_py = load('py_sat_sig10');

tmsv_fiber = fiber_py.data;
tmsv_sat = sat_py.data;

data_fiber = load('results_fiber_sig10-2');
data_sat = load('results_sat_sig10');

as_fiber = data_fiber.results;
as_sat = data_sat.results;



sq_fiber_py = load('py_sq_fiber_sig1-2');
sq_sat_py = load('py_sq_sat_sig1');

sq_tmsv_fiber = sq_fiber_py.data;
sq_tmsv_sat = sq_sat_py.data;

sq_data_fiber = load('sq_results_fiber_sig1-2');
sq_data_sat = load('sq_results_sat_sig1');

sq_as_fiber = sq_data_fiber.results;
sq_as_sat = sq_data_sat.results;


% Fiber
figure;
hold all;
set(0,'defaultTextInterpreter','latex');

Ls = linspace(30, 150, 24);

class_coh = 0.5*ones(length(Ls));
plot(Ls, class_coh, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

class_sq = 0.454*ones(length(Ls));
plot(Ls, class_sq, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

plot(Ls, tmsv_fiber, 'ko-','LineWidth', 1.2, 'DisplayName', 'TMSV $\left| \alpha \right\rangle$');
plot(Ls, as_fiber, 'v-', 'Color', [0 0.4470 0.7410],'LineWidth', 1.2, 'DisplayName', 'PA-PS $\left| \alpha \right\rangle$');

plot(Ls, sq_tmsv_fiber, 'ko--','LineWidth', 1.2, 'DisplayName', 'TMSV $\left| s \right\rangle$');
plot(Ls, sq_as_fiber, 'v--','Color', [0 0.4470 0.7410], 'LineWidth', 1.2, 'DisplayName', 'PA-PS $\left| s \right\rangle$');


hl = legend('show');
set(hl, 'Interpreter','latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$L$ [km]', 'Interpreter', 'latex');


txt1 = ('Classical limit  $\left| \alpha \right\rangle$');
txt2 = ('Classical limit  $\left| s \right\rangle$');
t1 = text(50.3, .490, txt1);
t2 = text(120, .444, txt2);

set(t1, 'Color','r');
set(t2, 'Color','r');

ylim([0.42 0.65]);
xlim([30 150]);

savefigures('fiber');

% Satellite
figure;
hold all;
set(gca, 'xscale', 'log')
set(0,'defaultTextInterpreter','latex');

zs = deg2rad([0 5 10 15 20 25 30 35 40 45 50 55 60, 65, 70]);


R = 6371000;
H = 500;
ds = sqrt((R*cos(zs)).^2 + H^2 + 2*R*H) - R * cos(zs);


% dmaxfit = 1250;
class_coh = 0.5*ones(200);
plot(linspace(500, max(ds), 200), class_coh, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

class_sq = 0.454*ones(200);
% plot(linspace(500, max(ds), 200), class_sq, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

plot(ds, tmsv_sat, 'ko-','LineWidth', 1.2, 'DisplayName', 'TMSV $\left| \alpha \right\rangle$');
plot(ds, as_sat, 'v-','Color', [0 0.4470 0.7410], 'LineWidth', 1.2, 'DisplayName', 'PA-PS $\left| \alpha \right\rangle$');


% plot(ds, sq_tmsv_sat, 'ko--','LineWidth', 1.2, 'DisplayName', 'TMSV $\left| \alpha \right\rangle$');
% plot(ds, sq_as_sat, 'o--','Color', [0 0.4470 0.7410], 'LineWidth', 1.2, 'DisplayName', 'PA-PS $\left| \alpha \right\rangle$');

% coeffs = polyfit(ds, as_sat, 4);
% Get fitted values
% fittedX = linspace(min(ds), dmaxfit, 200);
% fittedY = polyval(coeffs, fittedX);
% Plot the fitted line

% plot(fittedX, fittedY, 'b--', 'LineWidth', 0.5, 'DisplayName', 'PolyFit');

ylim([0.49 0.54])

hl = legend('show');
set(hl, 'Interpreter','latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$L$ [km]', 'Interpreter', 'latex');
savefigures('sat');

