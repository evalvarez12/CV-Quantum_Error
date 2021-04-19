clc;clear

OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');

% 'algorithm','sqp'
% e = 0.05;
me = (0.0086-.0033)/100;

tnum = 15;
% ts = linspace(0.01, 1, tnum);


% Fiber
Ls = linspace(50, 150, 20);

% Satellite
ts = [0.06143690562983698, 0.06086526299427021, 0.05917661503498977, 0.056443389085287135, 0.05323006348104734, 0.049082101315961746, 0.044150186248965786,  0.037171035021636226, 0.031114584577009226, 0.025515914750767258, 0.0198735248725103, 0.013815858165353754, 0.009422950626114226, 0.00582421, 0.00264689];
tnum = length(ts);

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


gmax = 10;
gmin = 0.1;
gini = 0.1;

% gmax = 1;
% gmin = 1;
% gini = 1;

Tmax = 0.999;
Tmin = 0.001;
Tini = 0.99;


rmax = 1;
rmin = 0;
rini = .5;

% sigma = [2, 5, 10, 20]
sig = 10;
eta = sqrt(10^(-1/10));
% eta = 0.7;

parfor i = 1:tnum
%     t = 10^(-tdB(i)/10);
%     L = Ls(i);
    disp([i]);
    
    % Fiber
%     t = 10^(-(0.16 * L)/10);
%     e = me * L + 0.0006;

    % Sat
    t = ts(i);
    e = t * 0.0186 + 0.0133/0.95;


%     fun_ep = @(par) -coh_loss_eta(par(1), par(2), t, e, 'epr', par(3), eta, sig);
%     fun_ps = @(par) -coh_loss_eta(par(1), par(2), t, e, 'ps', par(3), eta, sig);
%     fun_pa = @(par) -coh_loss_eta(par(1), par(2), t, e, 'pa', par(3), eta, sig);
%     fun_pc = @(par) -coh_loss_eta(par(1), par(2), t, e, 'pc', par(3), eta, sig);
    fun_as = @(par) -coh_loss_eta(par(1), par(2), t, e, 'as', par(3), eta, sig);
%     fun_sa = @(par) -coh_loss_eta(par(1), par(2), t, e, 'sa', par(3), eta, sig);
    
%     [~, Fepr(i)] = fmincon(fun_ep, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
%     [~, F_ps(i)] = fmincon(fun_ps, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
%     [~, F_pa(i)] = fmincon(fun_pa, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
%     [~, F_pc(i)] = fmincon(fun_pc, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
    [~, F_as(i)] = fmincon(fun_as, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
%     [~, F_sa(i)] = fmincon(fun_sa, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);


end

% 
% Fepr = -Fepr;
% F_ps = -F_ps;
% F_pa = -F_pa;
% F_pc = -F_pc;
F_as = -F_as; 
% F_sa = -F_sa;
% 
% Fepr
% F_ps
% F_pa
% F_pc
% F_as
% F_sa

% results = [Fepr(:), F_ps(:), F_pa(:), F_pc(:), F_as(:), F_sa(:)];
% save('results_sig20.mat', 'results');

results = F_as;

save('results_sat_sig10.mat', 'results');

