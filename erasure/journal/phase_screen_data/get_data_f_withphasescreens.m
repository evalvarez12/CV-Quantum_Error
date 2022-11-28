clear;
% Fidelity parameters
V = 3;
coh_sigma = 10;
% epsilon = 0.043;
eta = 0.9;



disp= [25.2491   25.2679   25.3129   25.4089   25.5921   25.9594   26.4423   27.1864   28.1080   29.2475   30.4720];
% disp= [25.4403   25.8218   26.6874   28.1420   30.0458   32.4353   34.8482   38.1391   41.3907   44.9236   48.2826];

erasure = true;

% Sample size - same as phase screen data size
N = 10000;

rec = '0.2';
dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000];

% rec = '0.15';
% dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600];

% rec = '0.1';
% dists = [1000 1200 1400 1600 1800 2000 2200 2400];

Fm = zeros(length(dists), N);
% Fstd = zeros(length(par), 1);


% Fm_dir = zeros(length(par), N);
% Fstd_dir = zeros(length(par), 1);

Fmax = zeros(N,1);
gmax = zeros(N, 1);
gs = -1.26:.01:1.26;

for i = 1:length(dists)

    epsilon = 10e-4 * disp(i) + 0.013;
    d = dists(i);

    T_ps = load(['data/ERASURE_d=', num2str(d), '_L0=1.5_l0=0.01_rec=', rec,'_10000.mat']);
    T_ps = T_ps.res;
    Ts1 = T_ps(randperm(length(T_ps)));
    Ts2 = T_ps(randperm(length(T_ps)))*logical(~erasure);
    Ts3 = T_ps(randperm(length(T_ps)));
    
    for j = 1:N
        T1 = Ts1(j);
        T2 = Ts2(j);
        T3 = Ts3(j);
    
        Fs = fid_tmsv_gen_loss_eq(V, gs, coh_sigma, T1, T2, T3, epsilon, eta);
        [Fmax(j), gmax(j)]= max(Fs);
    
    end
    
    
    Fm(i, :) = transpose(Fmax);
    Fstd(i) = std(Fmax);


    Fs_dir = fid_tmsv_dir_eq(Ts1, epsilon, coh_sigma);

    
    Fm_dir(i, :) = transpose(Fs_dir);
%     Fstd_dir(i) = std(Fs_dir);




end



if erasure
    save(['data/Fscatter_erasure_phasesc_rec=', rec, '_V', num2str(V),'.mat'], 'Fm');
    
else
    save(['data/Fscatter_phasesc_rec=', rec, '_V', num2str(V),'.mat'], 'Fm');

    save(['data/Fscatter_phasesc_rec=', rec, '_dir.mat'], 'Fm_dir');
end
