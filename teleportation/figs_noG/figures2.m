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

class = 0.5*ones(length(Ls));
plot(Ls, class, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

plot(Ls, tmsv_fiber, 'ko-','LineWidth', 1.2, 'DisplayName', 'TMSV');
plot(Ls, as_fiber, 'o-','LineWidth', 1.2, 'DisplayName', 'PA-PS');




legend();
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$L$ [km]', 'Interpreter', 'latex');

savefigures('fiber');

% Satellite
figure;
hold all;
set(gca, 'xscale', 'log')

zs = deg2rad([0 5 10 15 20 25 30 35 40 45 50 55 60, 65, 70]);


R = 6371000;
H = 500;
ds = sqrt((R*cos(zs)).^2 + H^2 + 2*R*H) - R * cos(zs);


% dmaxfit = 1250;
class = 0.5*ones(200);
plot(linspace(500, max(ds), 200), class, 'r-', 'LineWidth', 2, 'HandleVisibility','off');

plot(ds, tmsv_sat, 'ko-','LineWidth', 1.2, 'DisplayName', 'TMSV');
plot(ds, as_sat, 'o-','LineWidth', 1.2, 'DisplayName', 'PA-PS');


% coeffs = polyfit(ds, as_sat, 4);
% Get fitted values
% fittedX = linspace(min(ds), dmaxfit, 200);
% fittedY = polyval(coeffs, fittedX);
% Plot the fitted line

% plot(fittedX, fittedY, 'b--', 'LineWidth', 0.5, 'DisplayName', 'PolyFit');



legend();
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
xlabel('$L$ [km]', 'Interpreter', 'latex');
savefigures('sat');

