clear;
% Fidelity parameters
V = 3;
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
sigma = 0.; 


% Sample size
% N = 2e6;
N = 1;

par = linspace(0.4, 1, 24);

F_noerased = zeros(length(par), N);
F_1erased = zeros(length(par), N);
F_2erased = zeros(length(par), N);
F_3erased = zeros(length(par), N);

F_dir = zeros(length(par), N);

Fmax_noerased = zeros(N,1);
Fmax_1erased = zeros(N,1);
Fmax_2erased = zeros(N,1);
Fmax_3erased = zeros(N,1);

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
           
        Fs_noerased = fid_tmsv_gen_loss_eq(V, gs, coh_sigma, T1, T2, T3, epsilon);
        Fs_1erased = fid_tmsv_gen_loss_eq(V, gs, coh_sigma, 0, T2, T3, epsilon);
        Fs_2erased = fid_tmsv_gen_loss_eq(V, gs, coh_sigma, T1, 0, T3, epsilon);
        Fs_3erased = fid_tmsv_gen_loss_eq(V, gs, coh_sigma, T1, T2, 0, epsilon);
        
        [Fmax_noerased(j), ~]= max(Fs_noerased);
        [Fmax_1erased(j), ~]= max(Fs_1erased);
        [Fmax_2erased(j), ~]= max(Fs_2erased);
        [Fmax_3erased(j), ~]= max(Fs_3erased);
    end
    
    
    F_noerased(i, :) = transpose(Fmax_noerased);
    F_1erased(i, :) = transpose(Fmax_1erased);
    F_2erased(i, :) = transpose(Fmax_2erased);
    F_3erased(i, :) = transpose(Fmax_3erased);



    Fs_dir = fid_tmsv_dir_eq(Ts1, epsilon, coh_sigma);

    
    F_dir(i, :) = transpose(Fs_dir);


end



% mu - par = 0.4:0.002:1 sigma=0.2 pd=0 T0=0;
% sigma - par = 0:0.001:.3 mu=0.6 pd=0 T0=0;

results = cat(3, F_noerased,  F_1erased,  F_2erased, F_3erased, F_dir);


save('data/F_compare_V3_sigma10_sigma0.mat', 'results');
% save('data/F_pd_std.mat', 'Fstd');

% save('data/F_dir_sigma10.mat', 'Fm_dir');
% save('data/F_pd_std_dir.mat', 'Fstd_dir');
