clear;
% Fidelity parameters
V = 20;
coh_sigma = 10;
epsilon = 0;

% Distribution parameters
% Discrete error rate
pd = 0.;
% Discrete error transmissivity
T0 = 0.;
% Gaussian mean 
mu = 1;
% Gaussian sigma
sigma = 0.1; 


% Sample size
% N = 2e6;
N = 2e5;

par = linspace(0.4, 1, 60);

Fm = zeros(length(par), N);
% Fstd = zeros(length(par), 1);


% Fm_dir = zeros(length(par), N);
% Fstd_dir = zeros(length(par), 1);

Fmax = zeros(N,1);
gmax = zeros(N, 1);
gs = -1:.01:1;

for i = 1:length(par)

    var = par(i);

    pd = pd;
    T0 = T0;
    mu = var;
    sigma = sigma;

    Ts1 = distribution(pd, T0, mu, sigma, N);
    Ts2 = distribution(pd, T0, mu, sigma, N);
    Ts3 = distribution(pd, T0, mu, sigma, N);
    
    for j = 1:N
        T1 = Ts1(j);
        T2 = Ts2(j);
        T3 = Ts3(j);
    
        Fs = fid_tmsv_gen_loss_eq(V, gs, coh_sigma, T1, T2, T3, epsilon);
        [Fmax(j), gmax(j)]= max(Fs);
    
    end
    
    
    Fm(i, :) = transpose(Fmax);
%     Fstd(i) = std(Fmax);


%     Fs_dir = fid_tmsv_dir_eq(Ts2, epsilon, coh_sigma);

    
%     Fm_dir(i, :) = transpose(Fs_dir);
%     Fstd_dir(i) = std(Fs_dir);




end

% mu - par = 0.4:0.002:1 sigma=0.2 pd=0 T0=0;
% sigma - par = 0:0.001:.3 mu=0.6 pd=0 T0=0;

save('data/F_noerased_V20_sigma10.mat', 'Fm');
% save('data/F_pd_std.mat', 'Fstd');

% save('data/F_dir_sigma10.mat', 'Fm_dir');
% save('data/F_pd_std_dir.mat', 'Fstd_dir');
