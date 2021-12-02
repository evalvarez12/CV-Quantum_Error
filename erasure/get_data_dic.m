clc;clear

OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');


num = 24;
Ps = linspace(0, .7, num);

F = zeros(1, num);
Fdir = F;
pars = F;

sigma = 10;

Vmax = 7;
Vmin = 1;
Vini = 1.1;


dic = ParameterDic();
dic = dic.setupImpl();

F_123 = fid_tmsv(1, 1, 1, '123', sigma);
% F_123 = 0.0909 sigma = 10

parfor i = 1:num
    disp([i]);
    Pe = Ps(i);
    
    fun = @(par) -dic.Full_fidelity(Pe, par, sigma);

    [pars(i), F(i)] = fmincon(fun, Vini, [],[],[],[], Vmin, Vmax, [], OptOption);
        
    Fdir(i) = (1-Pe) + Pe * F_123;

    
end


F = -F;


results = [Ps(:), F(:), Fdir(:)];
save('data/results_F_dic.mat', 'results');
