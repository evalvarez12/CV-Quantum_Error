
V = 50;
sigma = 10;
epsilon = 0;

sigmaT = 0.;
meanT = 0.8;
T1r = normrnd(meanT, sigmaT);
T2r = normrnd(meanT, sigmaT);
T3r = normrnd(meanT, sigmaT);

% T1r = 0.9; 
% T2r = 0.9;
% T3r = 0.1;
disp(['T1 - ', num2str(T1r)]);
disp(['T2 - ', num2str(T2r)]);
disp(['T3 - ', num2str(T3r)]);

[minT, i] = min([T1r, T3r, T2r]);

g = i - 2;
g=0;
fid_tmsv_gen_loss_gs(V, g, sigma, T1r, T2r, T3r, epsilon)/2

T_mean = (T1r + T2r + T3r)/3;
%%%%% Direct
F_dir1 =  2 / (2*sigma*(1-sqrt(T1r))^2 + 2 + eps);
F_dir2 =  2 / (2*sigma*(1-sqrt(T2r))^2 + 2 + eps);
F_dir3 =  2 / (2*sigma*(1-sqrt(T3r))^2 + 2 + eps);

F_dir = (F_dir1 + F_dir2 + F_dir3)/3;
F_dir


exp(1).^((-1/4).*r.*((2+epsilon).*r+(sqrt(-1)*4).*2.^(1/2).*ri.*(( ...
  -1)+T1r.^(1/2)).*cos(theta)+(sqrt(-1)*(-4)).*2.^(1/2).*ra.*((-1)+ ...
  T1r.^(1/2)).*sin(theta))).*pi.^(-1).*r



% 
% % fid_tmsv_gen_loss(V, g, sigma, T1r, T2r, T3r, epsilon)
% fun1 = @(par) -fid_tmsv_gen_loss(V, par, sigma, T1r, T2r, T3r, epsilon);
% fun2 = @(par) -fid_tmsv_gen_loss(par, g, sigma, T1r, T2r, T3r, epsilon);
% 
% % fun3 = @(par) -fid_tmsv_gen_loss(par(1), par(2), sigma, T1r, T2r, T3r, epsilon);
% 
% pars=zeros(2:1);
% % F = 0;
% 
% gmax = 10;
% gmin = -10;
% gini = 1;
% 
% Vmax = 25;
% Vmin = 1;
% Vini = 1.5;
% 
% OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');
% 
% 
% % [pars(:,1), F] = fmincon(fun3, [gini, Vini], [],[],[],[], [gmin, Vmin], [gmax, Vmax], [], OptOption);
% 
% [g, F] = fmincon(fun1, gini, [],[],[],[], gmin, gmax, [], OptOption);
% % [V, F] = fmincon(fun2, Vini, [],[],[],[], Vmin, Vmax, [], OptOption);
% 
% 
% disp(['T1 - ', num2str(T1r)]);
% disp(['T2 - ', num2str(T2r)]);
% disp(['T3 - ', num2str(T3r)]);
% disp(['g - ', num2str(g)]);
% disp(['V - ', num2str(V)]);
% disp(['F - ', num2str(-F/2)]);
% disp('-----------------------');
