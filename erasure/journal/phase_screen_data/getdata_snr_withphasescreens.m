clear;
% Fidelity parameters
V = 10;
coh_sigma = 10;
epsilon = 0;



% Sample size - same as phase screen data size
N = 10000;

rec = '0.2';
dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000];

delta = 5;

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
    Ts1 = T_ps(randperm(length(T_ps)));
    Ts2 = T_ps(randperm(length(T_ps)));
    Ts3 = T_ps(randperm(length(T_ps)));
    
    for j = 1:N
        T1 = Ts1(j);
        T2 = Ts2(j);
        T3 = Ts3(j);
    
        snrs(j) = 10*log10(my_snr(T1, T2, T3, coh_sigma, V, delta));
    end

    snr_m(i) = mean(snrs);
    snr_std(i) = std(snrs);
    
end




save(['data/snrm_phasesc_rec=', rec, '_delta=', num2str(delta), '.mat'], 'snr_m');
save(['data/snrstd_phasesc_rec=', rec, '_delta=', num2str(delta), '.mat'], 'snr_std');


