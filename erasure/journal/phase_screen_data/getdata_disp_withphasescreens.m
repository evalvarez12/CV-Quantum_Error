clear;
% Fidelity parameters
V = 10;
coh_sigma = 10;
epsilon = 0;
eta = .9;
eps = 0.043;


% Sample size - same as phase screen data size
N = 10000;

rec = '0.1';
dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000];

% rec = '0.15';
% dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600];

% rec = '0.1';
% dists = [1000 1200 1400 1600 1800 2000 2200 2400];

snr_m = zeros(length(dists), 1);
snr_std = zeros(length(dists), 1);
snrs = zeros(N, 1);


for i = 1:length(dists)

    d = dists(i);

    T_ps = load(['data/ERASURE_d=', num2str(d), '_L0=1.5_l0=0.01_rec=', rec,'_10000.mat']);
    T_ps = T_ps.res;
    Ts1 = T_ps(randperm(length(T_ps)))*0;
    Ts2 = T_ps(randperm(length(T_ps)));
    Ts3 = T_ps(randperm(length(T_ps)));
    
    for j = 1:N
        T1 = Ts1(j);
        T2 = Ts2(j);
        T3 = Ts3(j);
    
        snrs(j) = bpsk_disp(T1, T2, T3, coh_sigma, V, eps, eta);
    end

    snr_m(i) = mean(snrs);
    snr_std(i) = std(snrs);
    
end




save(['data/dispm_phasesc_rec=', rec, '.mat'], 'snr_m');
save(['data/dispstd_phasesc_rec=', rec, '.mat'], 'snr_std');


