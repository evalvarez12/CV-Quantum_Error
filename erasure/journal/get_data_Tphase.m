clear;

% Generate PDF
fit_data = load('hist_data.mat');
fit_data = fit_data.c;

[ph_pdf, ~] = createFit(fit_data(:,1), fit_data(:,2));



pop = linspace(0, 1, 1e5);
w = ph_pdf(pop);
w(w<0)=0;
% Sample size
N= 1e5;


Ts1 = randsample(pop, N, true, w);
Ts2 = randsample(pop, N, true, w);
Ts3 = randsample(pop, N, true, w);


% Fidelity parameters
V = 20;
coh_sigma = 10;
epsilon = 0;


Fmax = zeros(N,1);
gmax = zeros(N, 1);
gs = -1:.01:1;

    
for j = 1:N
    T1 = Ts1(j);
    T2 = Ts2(j);
    T3 = Ts3(j);
    
    Fs = fid_tmsv_gen_loss_eq(V, gs, coh_sigma, T1, T2, T3, epsilon);
    [Fmax(j), gmax(j)]= max(Fs);
    
end
    

Fdir = fid_tmsv_dir_eq(Ts2, epsilon, coh_sigma);

figure;
hold on;
histogram(Fmax, 100);
histogram(Fdir, 100);


% mu - par = 0.4:0.002:1 sigma=0.2 pd=0 T0=0;
% sigma - par = 0:0.001:.3 mu=0.6 pd=0 T0=0;

% save('data/F_noerased_V20_sigma10.mat', 'Fm');
% save('data/F_pd_std.mat', 'Fstd');

% save('data/F_dir_sigma10.mat', 'Fm_dir');
% save('data/F_pd_std_dir.mat', 'Fstd_dir');
