clc;clear

OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');



tnum = 1;
F_ps = zeros(1, tnum);
F_pa = F_ps;
F_pc = F_ps;
F_as = F_ps;
F_sa = F_ps;
Fepr = F_ps;

opt_ps = F_ps;
opt_pa = F_ps;
opt_pc = F_ps;
opt_sa = F_ps;
opt_as = F_ps;

eta = sqrt(10^(-1/10));

gmax = 10;
gmin = eta;
gini = 1;

% gmax = 1;
% gmin = 1;
% gini = 1;

Tmax = 0.999;
Tmin = 0.001;
Tini = 0.99;


rmax = 1;
rmin = 0;
rini = .5;

sig = 0.575;
e = 0;
t = 0.6;

fun_ep = @(par) -sq_loss_eta(par(1), par(2), t, e, 'epr', par(3), eta, sig);
% fun_ps = @(par) -sq_loss_eta(par(1), par(2), t, e, 'ps', par(3), eta, sig);
% fun_pa = @(par) -sq_loss_eta(par(1), par(2), t, e, 'pa', par(3), eta, sig);
% fun_pc = @(par) -sq_loss_eta(par(1), par(2), t, e, 'pc', par(3), eta, sig);
% fun_as = @(par) -sq_loss_eta(par(1), par(2), t, e, 'as', par(3), eta, sig);
% fun_sa = @(par) -sq_loss_eta(par(1), par(2), t, e, 'sa', par(3), eta, sig);
    
    
    
[x, Fepr(1)] = fmincon(fun_ep, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
% [~, F_ps(1)] = fmincon(fun_ps, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
% [~, F_pa(1)] = fmincon(fun_pa, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
% [~, F_pc(1)] = fmincon(fun_pc, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
% [~, F_as(1)] = fmincon(fun_as, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
% [~, F_sa(1)] = fmincon(fun_sa, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);




Fepr = -Fepr;
% F_ps = -F_ps;
% F_pa = -F_pa;
% F_pc = -F_pc;
% F_as = -F_as; 
% F_sa = -F_sa;

Fepr
x
% F_ps
% F_pa
% F_pc
% F_as
% F_sa
r = 0.5407;
T = 0.069;
g = .1;
sq_loss_eta(r, T, t, e, 'epr', g, eta, sig)


