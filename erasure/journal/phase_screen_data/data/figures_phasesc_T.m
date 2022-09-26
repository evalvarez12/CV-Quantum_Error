clear;
set(0,'defaultTextInterpreter','latex');


% Fidelity parameters
V = 1;
coh_sigma = 10;
epsilon = 0;



% Sample size - same as phase screen data size
N = 10000;
rec = '0.1';
dists1 = [1000 1200 1400 1600 1800 2000 2200 2400 ];

Tm1 = zeros(length(dists1), 1);
Tstd1 = zeros(length(dists1), 1);

for i = 1:length(dists1)

    d = dists1(i);

    Ts = load(['data/ERASURE_d=', num2str(d), '_L0=1.5_l0=0.01_rec=', rec,'_10000.mat']);
    Ts = Ts.res;

    Tm1(i) = mean(Ts);
    Tstd1(i) = std(Ts);

end


rec = '0.15';
dists2 = [1000 1200 1400 1600 1800 2000 2200 2400 2600];


Tm2 = zeros(length(dists2), 1);
Tstd2 = zeros(length(dists2), 1);
for i = 1:length(dists2)

    d = dists2(i);

    Ts = load(['data/ERASURE_d=', num2str(d), '_L0=1.5_l0=0.01_rec=', rec,'_10000.mat']);
    Ts = Ts.res;

    Tm2(i) = mean(Ts);
    Tstd2(i) = std(Ts);
end




Tm1 = transpose(Tm1);
Tstd1 = transpose(Tstd1);


y1 = [Tm1(1:end-1) - Tstd1(1:end-1); Tm1(2:end) - Tstd1(2:end); flipud(Tm1(2:end) + Tstd1(2:end)); flipud(Tm1(1:end-1) + Tstd1(1:end-1))];

x1 = [dists1(1:end-1); dists1(2:end); dists1(2:end); dists1(1:end-1)];


Tm2 = transpose(Tm2);
Tstd2 = transpose(Tstd2);

y2 = [Tm2(1:end-1) - Tstd2(1:end-1); Tm2(2:end) - Tstd2(2:end); flipud(Tm2(2:end) + Tstd2(2:end)); flipud(Tm2(1:end-1) + Tstd2(1:end-1))];

x2 = [dists2(1:end-1); dists2(2:end); dists2(2:end); dists2(1:end-1)];

hold on;
fill(x1, y1, 'blue','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');

plot(dists1, Tm1, 'DisplayName', '$r_d=0.1m$');


fill(x2, y2, 'red','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');

plot(dists2, Tm2, 'DisplayName', '$r_d=0.15m$');

% errorbar(par(1:2:end), snr_m(1:2:end), snr_std(1:2:end),'-','LineWidth',1.2, 'CapSize',4, 'DisplayName','\Delta =10');



% errorbar(dists, Tm, Tstd, 'DisplayName', 'r_d=0.15m')





legend('Interpreter', 'latex')
xlabel('$L[m]$', 'Interpreter', 'latex');
ylabel('$T$', 'Interpreter', 'latex');
savefigures('phasesc_T')
