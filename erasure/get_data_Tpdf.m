clear;
% Fidelity parameters
V = 20;
coh_sigma = 10;
epsilon = 0;

% Distribution parameters
% Discrete error rate
pd = 0.1;
% Discrete error transmissivity
T0 = 0.;
% Gaussian mean 
mu = 1;
% Gaussian sigma
sigma = 0.0; 


% Sample size
% N = 2e6;
N = 2e3;

par = 0:0.005:.5;

Fm = zeros(length(par), 1);
Fstd = zeros(length(par), 1);


Fm_dir = zeros(length(par), 1);
Fstd_dir = zeros(length(par), 1);

Fmax = zeros(N,1);
gmax = zeros(N, 1);
gs = -1:.01:1;

for i = 1:length(par)

    var = par(i);

    Ts1 = distribution(var, T0, mu, sigma, N);
    Ts2 = distribution(var, T0, mu, sigma, N);
    Ts3 = distribution(var, T0, mu, sigma, N);
    
    for j = 1:N
        T1 = Ts1(j);
        T2 = Ts2(j);
        T3 = Ts3(j);
    
        Fs = fid_tmsv_gen_loss_eq(V, gs, coh_sigma, T1, T2, T3, epsilon);
        [Fmax(j), gmax(j)]= max(Fs);
    
    end
    
    
    Fm(i) = mean(Fmax);
    Fstd(i) = std(Fmax);


    Fs_dir = fid_tmsv_dir_eq(Ts1, epsilon, coh_sigma);

    
    Fm_dir(i) = mean(Fs_dir);
    Fstd_dir(i) = std(Fs_dir);




end

save('data/F_pd_mean.mat', 'Fm');
save('data/F_pd_std.mat', 'Fstd');

save('data/F_pd_mean_dir.mat', 'Fm_dir');
save('data/F_pd_std_dir.mat', 'Fstd_dir');
