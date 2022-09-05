V = 10;
sigma = 10;

g = -1;

Ts = linspace(0.4,1,12);
l = length(Ts);

F_tmsv = zeros(1, l);
F_ps = zeros(1, l);
F_pa = zeros(1, l);
F_pc = zeros(1, l);
F_pspa = zeros(1, l);
F_paps = zeros(1, l);

Tmax = 0.999;
Tmin = 0.001;
Tini = 0.99;

OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');


for i=1:l
    
    T = Ts(i);
    T1 = 0;
    T2 = T;
    T3 = T;

    fun_ps = @(par) -f_code(T1, T2, T3, V, par, g, 'ps', sigma);
    fun_pa = @(par) -f_code(T1, T2, T3, V, par, g, 'pa', sigma);
    fun_pc = @(par) -f_code(T1, T2, T3, V, par, g, 'pc', sigma);
    fun_pspa = @(par) -f_code(T1, T2, T3, V, par, g, 'pspa', sigma);
    fun_paps = @(par) -f_code(T1, T2, T3, V, par, g, 'paps', sigma);

    F_tmsv(i) = f_code(T1, T2, T3, V, 0, g, 'tmsv', sigma);
    [~, F_ps(i)] = fmincon(fun_ps, Tini, [],[],[],[], Tmin, Tmax, [], OptOption);
    [~, F_pa(i)] = fmincon(fun_pa, Tini, [],[],[],[], Tmin, Tmax, [], OptOption);
    [~, F_pc(i)] = fmincon(fun_pc, Tini, [],[],[],[], Tmin, Tmax, [], OptOption);
    [~, F_pspa(i)] = fmincon(fun_pspa, Tini, [],[],[],[], Tmin, Tmax, [], OptOption);
    [~, F_paps(i)] = fmincon(fun_paps, Tini, [],[],[],[], Tmin, Tmax, [], OptOption);
end

F_ps = -F_ps;
F_pa = -F_pa;
F_pc = -F_pc;
F_pspa = -F_pspa;
F_paps = -F_paps;

results = [Ts(:), F_tmsv(:), F_ps(:), F_pa(:), F_pc(:), F_pspa(:), F_paps(:)];
save('fid_T1erased_V10_s10.mat', 'results');