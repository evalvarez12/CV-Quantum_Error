clc;clear

OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');

% 'algorithm','sqp'


num = 30;

F_as = zeros(num, num);
ts = linspace(0.0001, 1, num);
es = logspace(-2.2, -0.9, 30);


gmax = 10;
gmin = 0.1;
gini = 0.1;


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

parfor i = 1:num
    for j = 1:num
        disp(['---->', num2str(i), ' ', num2str(j)]);
        t = ts(i);
        e = es(j);
       
        fun_as = @(par) -coh_loss_eta(par(1), par(2), t, e, 'as', par(3), eta, sig);
        
        [~, F_as(i, j)] = fmincon(fun_as, [Tini, rini, gini], [],[],[],[], [Tmin, rmin, gmin], [Tmax, rmax, gmax], [], OptOption);

    end
F_
end


F_as = -F_as; 



results = F_as;

save('data-as-sig10.mat', 'results');

