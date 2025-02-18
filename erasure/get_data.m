clc;clear

OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');


V_num = 72;
Vs = linspace(1, 9, V_num);
% Vs = [3];

F_1 = zeros(1, V_num);
F_2 = F_1;
F_3 = F_1;
F_12 = F_1;
F_13 = F_1;
F_23 = F_1;

par1 = zeros(1, V_num);
par2 = par1;
par3 = par1;
par12 = par1;
par13 = par1;
par23 = par1;

sigma = 10;

gmax = 10;
gmin = -10;
gini = 1;

% gmax = 1;
% gmin = 1;
% gini = 1;

Tmax = 0.999;
Tmin = 0.001;
Tini = 0.99;


Vmax = 10;
Vmin = 1;
Vini = 1.5;




% fun_1 = @(par) -fid_tmsv(par(1), par(2), par(3), '1', sigma);
% fun_2 = @(par) -fid_tmsv(par(1), par(2), par(3), '2', sigma);
% fun_3 = @(par) -fid_tmsv(par(1), par(2), par(3), '3', sigma);
% fun_12 = @(par) -fid_tmsv(par(1), par(2), par(3), '12', sigma);
% fun_13 = @(par) -fid_tmsv(par(1), par(2), par(3), '13', sigma);
% fun_23 = @(par) -fid_tmsv(par(1), par(2), par(3), '23', sigma);
% 
% 
% [~, F_1] = fmincon(fun_ep, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
% [~, F_2] = fmincon(fun_ps, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
% [~, F_3] = fmincon(fun_pa, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
% [~, F_12] = fmincon(fun_pc, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
% [~, F_13] = fmincon(fun_as, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);
% [~, F_23 ] = fmincon(fun_sa, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);


for i = 1:V_num
    V = Vs(i);
%     
    fun_1 = @(par) -fid_tmsv_ri(V, par(1), '1', sigma);
    fun_2 = @(par) -fid_tmsv_ri(V, par(1), '2', sigma);
    fun_3 = @(par) -fid_tmsv_ri(V, par(1), '3', sigma);
    fun_12 = @(par) -fid_tmsv_ri(V, par(1), '12', sigma);
    fun_13 = @(par) -fid_tmsv_ri(V, par(1), '13', sigma);
    fun_23 = @(par) -fid_tmsv_ri(V, par(1), '23', sigma);


    [par1(i), F_1(i)] = fmincon(fun_1, gini, [],[],[],[], gmin, gmax, [], OptOption);
    [par2(i), F_2(i)] = fmincon(fun_2, gini, [],[],[],[], gmin, gmax, [], OptOption);
    [par3(i), F_3(i)] = fmincon(fun_3, gini, [],[],[],[], gmin, gmax, [], OptOption);
    [par12(i), F_12(i)] = fmincon(fun_12, gini, [],[],[],[], gmin, gmax, [], OptOption);
    [par13(i), F_13(i)] = fmincon(fun_13, gini, [],[],[],[], gmin, gmax, [], OptOption);
    [par23(i), F_23(i)] = fmincon(fun_23, gini, [],[],[],[], gmin, gmax, [], OptOption);
    

end


F_1 = -F_1;
F_2 = -F_2;
F_3 = -F_3;
F_12 = -F_12;
F_13 = -F_13; 
F_23 = -F_23;



results = [Vs(:), F_1(:), F_2(:), F_3(:), F_12(:), F_13(:), F_23(:)];
save('data/dic_Fs.mat', 'results');

results_pars = [transpose(par1(:)) ...
                transpose(par2(:)) ...
                transpose(par3(:)) ...
                transpose(par12(:)) ...
                transpose(par13(:)) ...
                transpose(par23(:))];
save('data/dic_Fs_pars.mat', 'results_pars');

