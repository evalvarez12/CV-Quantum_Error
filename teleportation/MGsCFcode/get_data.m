clc;clear

OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');

% 'algorithm','sqp'
% e = 0.05;
me = (0.0086-.0033)/100;

% tnum = 15;
% ts = linspace(0.01, 1, tnum);


% Fiber
% Ls = linspace(50, 150, 20);

% Satellite
ts = [0.06143690562983698, 0.06086526299427021, 0.05917661503498977, 0.056443389085287135, 0.05323006348104734, 0.049082101315961746, 0.044150186248965786,  0.037171035021636226, 0.031114584577009226, 0.025515914750767258, 0.0198735248725103, 0.013815858165353754, 0.009422950626114226, 0.00582421, 0.00264689];
tnum = length(ts);
zs = deg2rad([0 5 10 15 20 25 30 35 40 45 50 55 60 65 70]);


R = 6371000;
H = 500;
ds = sqrt((R*cos(zs)).^2 + H^2 + 2*R*H) - R * cos(zs);



F_ps = zeros(1, tnum);
F_pa = F_ps;
F_pc = F_ps;
F_as = F_ps;
F_sa = F_ps;
Fepr = F_ps;

F_as1 = F_ps;
F_as2 = F_ps;

eta = sqrt(10^(-1/10));
% s = 0.575;
sig = 10;

gmax = 10;
gmin = 0.1;
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

r1 = 0.230;
r2 = 0.575;

s = 0.25;


for i = 6:tnum
    disp([i]);
    
    % Fixed
%     t = ts(i);
%     e = 0.05


    % Fiber
%     t = 10^(-tdB(i)/10);
%     L = Ls(i);
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
    

%     fun_as1 = @(par) -coh_loss_eta(par(1), r1, t, e, 'as', 1/eta, eta, sig);
%     fun_as2 = @(par) -coh_loss_eta(par(1), r2, t, e, 'as', 1/eta, eta, sig);

    
%     fun_ep = @(par) -sq_loss_eta(par(1), par(2), t, e, 'epr', par(3), eta, sig);
%     fun_ps = @(par) -sq_loss_eta(par(1), par(2), t, e, 'ps', par(3), eta, sig);
%     fun_pa = @(par) -sq_loss_eta(par(1), par(2), t, e, 'pa', par(3), eta, sig);
%     fun_pc = @(par) -sq_loss_eta(par(1), par(2), t, e, 'pc', par(3), eta, sig);
%     fun_as = @(par) -sq_loss_eta(par(1), par(2), t, e, 'as', par(3), eta, sig);
%     fun_sa = @(par) -sq_loss_eta(par(1), par(2), t, e, 'sa', par(3), eta, sig);
    
%     fun_ep = @(par) -sq_coh_loss_eta(par(1), par(2), t, e, 'epr', par(3), eta, s, sig);
%     fun_ps = @(par) -sq_coh_loss_eta(par(1), par(2), t, e, 'ps', par(3), eta, s, sig);
%     fun_pa = @(par) -sq_coh_loss_eta(par(1), par(2), t, e, 'pa', par(3), eta, s, sig);
%     fun_pc = @(par) -sq_coh_loss_eta(par(1), par(2), t, e, 'pc', par(3), eta, s, sig);
%     fun_as = @(par) -sq_coh_loss_eta(par(1), par(2), t, e, 'as', par(3), eta, s, sig);
%     fun_sa = @(par) -sq_coh_loss_eta(par(1), par(2), t, e, 'sa', par(3), eta, s, sig);
%     
%     
%     [~, Fepr(i)] = fmincon(fun_ep, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
%     [~, F_ps(i)] = fmincon(fun_ps, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
%     [~, F_pa(i)] = fmincon(fun_pa, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
%     [~, F_pc(i)] = fmincon(fun_pc, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
    [~, F_as(i)] = fmincon(fun_as, [Tini], [],[],[],[], [Tmin], [Tmax], [], OptOption);
%     [~, F_sa(i)] = fmincon(fun_sa, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);


%     [~, F_as1(i)] = fmincon(fun_as1, [Tini], [],[],[],[], [Tmin], [Tmax], [], OptOption);
%     [~, F_as2(i)] = fmincon(fun_as2, [Tini], [],[],[],[], [Tmin], [Tmax], [], OptOption);

    
    
%     disp([ds(i), F_as(i)]);

end


% Fepr = -Fepr;
% F_ps = -F_ps;
% F_pa = -F_pa;
% F_pc = -F_pc;
F_as = -F_as; 
% F_sa = -F_sa;

% Fepr
% F_ps
% F_pa
% F_pc
F_as
% F_sa

% results = [Fepr(:), F_ps(:), F_pa(:), F_pc(:), F_as(:), F_sa(:)];
% save('sq_results_r1.mat', 'results');

results = [F_as(:), F_as1(:), F_as2(:)];

% save('results_sat_sig10.mat', 'results');

