clear;
% Fidelity parameters
V = 1;
coh_sigma = 10;
epsilon = 0;



% Sample size - same as phase screen data size
N = 10000;
rec = '0.15';
dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600];

Tm = zeros(length(dists), 1);
Tstd = zeros(length(dists), 1);



for i = 1:length(dists)

    d = dists(i);

    Ts = load(['data/ERASURE_d=', num2str(d), '_L0=1.5_l0=0.01_rec=', rec,'_10000.mat']);
    Ts = Ts.res;

    Tm(i) = mean(Ts);
    Tstd(i) = std(Ts);

end

errorbar(dists, Tm, Tstd, 'DisplayName', 'r_d=0.15m')

legend()
xlabel('$L[m]$', 'Interpreter', 'latex');
ylabel('$T$', 'Interpreter', 'latex');
savefigures('phasesc_T')
