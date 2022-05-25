clc;clear


V = 10;
g = -1;
sigma = 10;
T1r = 0.9;
T2r = 0.6;
T3r = 0.3;
epsilon = 0;

% fid_tmsv_gen_loss(V, g, sigma, T1r, T2r, T3r, epsilon)
fun1 = @(par) -fid_tmsv_gen_loss(V, par, sigma, T1r, T2r, T3r, epsilon);
fun2 = @(par) -fid_tmsv_gen_loss(par, g, sigma, T1r, T2r, T3r, epsilon);

% fun3 = @(par) -fid_tmsv_gen_loss(par(1), par(2), sigma, T1r, T2r, T3r, epsilon);

pars=zeros(2:1);
% F = 0;

gmax = 10;
gmin = -10;
gini = 1;

Vmax = 25;
Vmin = 1;
Vini = 1.5;

OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');

% [pars(:,1), F] = fmincon(fun3, [gini, Vini], [],[],[],[], [gmin, Vmin], [gmax, Vmax], [], OptOption);

[g_opt, F] = fmincon(fun1, gini, [],[],[],[], gmin, gmax, [], OptOption);
% [g_opt, F] = fmincon(fun2, Vini, [],[],[],[], Vmin, Vmax, [], OptOption);

-F/2
g_opt