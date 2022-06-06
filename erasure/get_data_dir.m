
% Parameters
sigma = 10;
epsilon = 0;

% Sample size
N = 2000000;

% Distribution properties
sigmaT = 0.3;
meanT = 0.6;


Ts = normrnd(meanT, sigmaT, N, 1); 
Ts = max(Ts, 0);
Ts(Ts>1) = 1;
Fs = fid_tmsv_dir_eq(Ts, epsilon, sigma);

mean(Fs)


% T = 0.6 s = 0.1 - 0.6633
% T = 0.6 s = 0.2 - 0.6578
% T = 0.6 s = 0.3 - 0.6442