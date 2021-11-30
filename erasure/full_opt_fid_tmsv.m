clc;clear

OptOption = optimoptions('fmincon', 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');

sigma = 10;

gmax = 10;
gmin = -10; 
gini = 1;

Vmax = 100;
Vmin = 1;
Vini = 1.5;


fun_1 = @(par) -fid_tmsv(par(1), par(2), par(3), '1', sigma);
fun_2 = @(par) -fid_tmsv(par(1), par(2), par(3), '2', sigma);
fun_3 = @(par) -fid_tmsv(par(1), par(2), par(3), '3', sigma);
fun_12 = @(par) -fid_tmsv(par(1), par(2), par(3), '12', sigma);
fun_13 = @(par) -fid_tmsv(par(1), par(2), par(3), '13', sigma);
fun_23 = @(par) -fid_tmsv(par(1), par(2), par(3), '23', sigma);


[~, F_1] = fmincon(fun_1, [Vini, gini, gini], [],[],[],[], [Vmin, gmin, gmin], [Vmax, gmax, gmax], [], OptOption);
[~, F_2] = fmincon(fun_2, [Vini, gini, gini], [],[],[],[], [Vmin, gmin, gmin], [Vmax, gmax, gmax], [], OptOption);
[~, F_3] = fmincon(fun_3, [Vini, gini, gini], [],[],[],[], [Vmin, gmin, gmin], [Vmax, gmax, gmax], [], OptOption);
[~, F_12] = fmincon(fun_12, [Vini, gini, gini], [],[],[],[], [Vmin, gmin, gmin], [Vmax, gmax, gmax], [], OptOption);
[~, F_13] = fmincon(fun_13, [Vini, gini, gini], [],[],[],[], [Vmin, gmin, gmin], [Vmax, gmax, gmax], [], OptOption);
[~, F_23] = fmincon(fun_23, [Vini, gini, gini], [],[],[],[], [Vmin, gmin, gmin], [Vmax, gmax, gmax], [], OptOption);
F123 = fid_tmsv(1, 1, 1, '123', sigma);

F_1 = -F_1;
F_2 = -F_2;
F_3 = -F_3;
F_12 = -F_12;
F_13 = -F_13;
F_23 = -F_23;

num = 12;
Ps = linspace(0, .6, num);

F = zeros(1, num);
Fdir = F;


sigma = 10;
    

for i = 1:num
    disp([i]);
    Pe = Ps(i);
    
    
    F(i) = (1-Pe)^3  + Pe*(1-Pe)^2*F_1 + + Pe*(1-Pe)^2*F_2 + + Pe*(1-Pe)^2*F_3 + ...
    Pe^2*(1-Pe)*F_12 + Pe^2*(1-Pe)*F_13 +Pe^2*(1-Pe)*F_23 + Pe^3*F123;
    
    Fdir(i) = (1-Pe) + Pe * F123;
    
end

results = [F(:), Fdir(:)];
save('results_fullopt_Ftmsv.mat', 'results');