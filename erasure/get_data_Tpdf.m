
V = 50;
sigma = 10;
epsilon = 0;

sigmaT = 0.2;
meanT = 0.6;
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

disp([ 0, fid_tmsv_gen_loss_gs(V, 0, sigma, T1r, T2r, T3r, epsilon)/2]);
disp([ 1, fid_tmsv_gen_loss_gs(V, 1, sigma, T1r, T2r, T3r, epsilon)/2]);
disp([ -1, fid_tmsv_gen_loss_gs(V, -1, sigma, T1r, T2r, T3r, epsilon)/2]);

disp([ 0.5, fid_tmsv_gen_loss(V, 0.5, sigma, T1r, T2r, T3r, epsilon)/2]);
disp([ -0.5, fid_tmsv_gen_loss(V, -0.5, sigma, T1r, T2r, T3r, epsilon)/2]);
disp([ 1, fid_tmsv_gen_loss(V, 1, sigma, T1r, T2r, T3r, epsilon)/2]);
disp([ -1, fid_tmsv_gen_loss(V, -1, sigma, T1r, T2r, T3r, epsilon)/2]);



T_mean = (T1r + T2r + T3r)/3;
%%%%% Direct
F_dir1 =  2 / (2*sigma*(1-sqrt(T1r))^2 + 2 + eps);
F_dir2 =  2 / (2*sigma*(1-sqrt(T2r))^2 + 2 + eps);
F_dir3 =  2 / (2*sigma*(1-sqrt(T3r))^2 + 2 + eps);

F_dir = (F_dir1 + F_dir2 + F_dir3)/3;
F_dir




% T = 0.8;
% F_dir =  2 / (2*sigma*(1-sqrt(T))^2 + 2 + eps);
% F_dir*2
% 
% fid_tmsv_dir(T, 0, sigma)


