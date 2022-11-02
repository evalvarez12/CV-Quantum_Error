V = 10;
sigma = 10;

g = -1;

Ts = linspace(0.6,1,12);
len = length(Ts);

F_tmsv = zeros(1, len);
F_ps = zeros(1, len);
F_pa = zeros(1, len);
F_pc = zeros(1, len);
F_pspa = zeros(1, len);
F_paps = zeros(1, len);
F_sb = zeros(1,len);

Tmax = 0.999;
Tmin = 0.001;
Tini = 0.99;

dmax = pi;
dmin = 0.0;
dini = 0.4;

OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');


for i=11:len
    
    T = Ts(i);
    T1 = 0;
    T2 = T;
    T3 = T;

%     fun_ps = @(par) -f_code(T1, T2, T3, V, par, g, 'ps', sigma);
%     fun_pa = @(par) -f_code(T1, T2, T3, V, par, g, 'pa', sigma);
%     fun_pc = @(par) -f_code(T1, T2, T3, V, par, g, 'pc', sigma);
%     fun_pspa = @(par) -f_code(T1, T2, T3, V, par, g, 'pspa', sigma);
%     fun_paps = @(par) -f_code(T1, T2, T3, V, par, g, 'paps', sigma);
    fun_sb = @(par) -f_code_sb(T1, T2, T3, V, par, g, sigma);

    F_tmsv(i) = f_code(T1, T2, T3, V, 0, g, 'tmsv', sigma);
%     [~, F_ps(i)] = fmincon(fun_ps, Tini, [],[],[],[], Tmin, Tmax, [], OptOption);
%     [~, F_pa(i)] = fmincon(fun_pa, Tini, [],[],[],[], Tmin, Tmax, [], OptOption);
%     [~, F_pc(i)] = fmincon(fun_pc, Tini, [],[],[],[], Tmin, Tmax, [], OptOption);
%     [~, F_pspa(i)] = fmincon(fun_pspa, Tini, [],[],[],[], Tmin, Tmax, [], OptOption);
%     [~, F_paps(i)] = fmincon(fun_paps, Tini, [],[],[],[], Tmin, Tmax, [], OptOption);
    [d, F_sb(i)] = fmincon(fun_sb, dini, [],[],[],[], dmin, dmax, [], OptOption);
    d

end

F_ps = -F_ps;
F_pa = -F_pa;
F_pc = -F_pc;
F_pspa = -F_pspa;
F_paps = -F_paps;
F_sb = -F_sb;

results = [Ts(:), F_tmsv(:), F_ps(:), F_pa(:), F_pc(:), F_pspa(:), F_paps(:),F_sb(:)];
% save('fid_T1erased_V10_s10.mat', 'results');