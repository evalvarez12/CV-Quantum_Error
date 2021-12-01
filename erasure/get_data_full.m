clc;clear

OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');

% 'algorithm','sqp'
me = (0.0086-.0033)/100;


num = 12;
Ps = linspace(0, .6, num);

F = zeros(1, num);
pars = zeros(3, num);
Fdir = F;

sigma = 10;

gmax = 10;
gmin = -10;
gini = 1;



Vmax = 100;
Vmin = 1;
Vini = 1.5;

F123 = fid_tmsv(1, 1, 1, '123', sigma);



parfor i = 1:num
    disp([i]);
    Pe = Ps(i);
    
    fun = @(par) -full_fid_tmsv(par(1), par(2), par(3), Pe, sigma);

    [pars(:, i), F(i)] = fmincon(fun, [Vini, gini, gini], [],[],[],[], [Vmin, gmin, gmin], [Vmax, gmax, gmax], [], OptOption);
        
    Fdir(i) = (1-Pe) + Pe * F123;

    
end


F = -F;


results = [Ps(:), transpose(pars(1,:)), F(:), Fdir(:)];
save('results_Ftmsv_full.mat', 'results');
