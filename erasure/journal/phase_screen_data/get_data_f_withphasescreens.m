clear;
% Fidelity parameters
V = 3;
coh_sigma = 10;
epsilon = 0;



% Sample size - same as phase screen data size
N = 10000;

dists = [900 1100 1300 1500 1700 1900 2100 2300];

Fm = zeros(length(dists), N);
% Fstd = zeros(length(par), 1);


% Fm_dir = zeros(length(par), N);
% Fstd_dir = zeros(length(par), 1);

Fmax = zeros(N,1);
gmax = zeros(N, 1);
gs = -1:.01:1;

for i = 1:length(dists)

    d = dists(i);

    T_ps = load(['ERASURE_d=', num2str(d), '_L0=1.5_l0=0.01_10000.mat']);
    T_ps = T_ps.res;
    Ts1 = T_ps(randperm(length(T_ps)));
    Ts2 = T_ps(randperm(length(T_ps)));
    Ts3 = T_ps(randperm(length(T_ps)));
    
    for j = 1:N
        T1 = Ts1(j);
        T2 = Ts2(j);
        T3 = Ts3(j);
    
        Fs = fid_tmsv_gen_loss_eq(V, gs, coh_sigma, T1, T2, T3, epsilon);
        [Fmax(j), gmax(j)]= max(Fs);
    
    end
    
    
    Fm(i, :) = transpose(Fmax);
    Fstd(i) = std(Fmax);


    Fs_dir = fid_tmsv_dir_eq(Ts1, epsilon, coh_sigma);

    
    Fm_dir(i, :) = transpose(Fs_dir);
%     Fstd_dir(i) = std(Fs_dir);




end

% mu - par = 0.4:0.002:1 sigma=0.2 pd=0 T0=0;
% sigma - par = 0:0.001:.3 mu=0.6 pd=0 T0=0;

save('data/Fscatter_phasesc_V3.mat', 'Fm');
% save('data/F_pd_std.mat', 'Fstd');

save('data/Fscatter_phasesc_dir.mat', 'Fm_dir');
% save('data/F_pd_std_dir.mat', 'Fstd_dir');
