
% Parameters
sigma = 10;
epsilon = 0;

% Sample size
N = 2e6;

meanT = 0.6;
sigmaT = 0:0.01:0.5;

betaT = 1:0.03:3;

Fm = zeros(length(betaT), 1);
Fstd = zeros(length(betaT), 1);

for i = 1:length(betaT)


%     Ts = normrnd(meanT, sigmaT(i), N, 1); 
    Ts = betarnd(4, betaT(i), N, 1);

    Ts = max(Ts, 0);
    Ts(Ts>1) = 1;
    Fs = fid_tmsv_dir_eq(Ts, epsilon, sigma);

    
    Fm(i) = mean(Fs);
    Fstd(i) = std(Fs);

end

% meanT = 0.6;
% sigmaT = 0:0.01:0.5;
% save('data/F_mean_dir.mat', 'Fm');
% save('data/F_std_dir.mat', 'Fstd');

% betaT = 1:0.03:3; ---- 4
save('data/F_beta_mean_dir.mat', 'Fm');
save('data/F_beta_std_dir.mat', 'Fstd');
