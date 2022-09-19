clear all;
close all;
set(0,'defaultTextInterpreter','latex');

% Parameters
V = 20;
coh_sigma = 10;
epsilon = 0;
delta = 10;

% Distribution parameters
% Discrete error rate
pd = 0.;
% Discrete error transmissivity
T0 = 0.;
% Gaussian mean 
mu = 0.6;
% Gaussian sigma
sigma = 0.2; 


% Sample size
% N = 2e6;
N = 2e6;

par = 0.4:0.01:1;

snr_m = zeros(length(par), 1);
snr_std = zeros(length(par), 1);
snr = zeros(N, 1);


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
    
       snr(j) = SNR(T1, T2, T3, sigma, V, delta);

    
    end
    
    
   snr_m(i) = mean(snr);
   snr_std(i) = std(snr);
end

save('data/snr_m_t3.mat', 'snr_m');
save('data/snr_std_t3.mat', 'snr_std');



