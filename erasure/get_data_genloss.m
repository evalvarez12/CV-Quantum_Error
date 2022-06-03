

V = 1;
g = -1;
sigma = 10;
T1r = 0.7;
T2r = 0.7;
T3r = 0.7;
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

[g, F] = fmincon(fun1, gini, [],[],[],[], gmin, gmax, [], OptOption);
% [V, F] = fmincon(fun2, Vini, [],[],[],[], Vmin, Vmax, [], OptOption);


disp(['T1 - ', num2str(T1r)]);
disp(['T2 - ', num2str(T2r)]);
disp(['T3 - ', num2str(T3r)]);
disp(['g - ', num2str(g)]);
disp(['V - ', num2str(V)]);
disp(['F - ', num2str(-F/2)]);
disp('-----------------------');
