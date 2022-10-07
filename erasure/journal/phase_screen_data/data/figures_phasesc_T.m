clear;
set(0,'defaultTextInterpreter','latex');


% Sample size - same as phase screen data size
N = 10000;
rec = '0.1';
dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000];


Tm1 = zeros(length(dists), 1);
Tstd1 = zeros(length(dists), 1);
sigs1 = zeros(length(dists), 1);

for i = 1:length(dists)

    d = dists(i);

    Ts = load(['data/ERASURE_d=', num2str(d), '_L0=1.5_l0=0.01_rec=', rec,'_10000.mat']);
    Ts = Ts.res;

    Tm1(i) = mean(Ts);
    Tstd1(i) = std(Ts);

    sigs1(i) = (mean(Ts.^2)/(mean(Ts)^2)) - 1;

end


rec = '0.15';


Tm2 = zeros(length(dists), 1);
Tstd2 = zeros(length(dists), 1);
sigs2 = zeros(length(dists), 1);

for i = 1:length(dists)

    d = dists(i);

    Ts = load(['data/ERASURE_d=', num2str(d), '_L0=1.5_l0=0.01_rec=', rec,'_10000.mat']);
    Ts = Ts.res;
%     Ts = -10*log10(Ts.res);

    Tm2(i) = mean(Ts);
    Tstd2(i) = std(Ts);

    sigs2(i) = (mean(Ts.^2)/(mean(Ts)^2)) - 1;

end

rec = '0.2';


Tm3 = zeros(length(dists), 1);
Tstd3 = zeros(length(dists), 1);
sigs3 = zeros(length(dists), 1);

for i = 1:length(dists)

    d = dists(i);

    Ts = load(['data/ERASURE_d=', num2str(d), '_L0=1.5_l0=0.01_rec=', rec,'_10000.mat']);
    Ts = Ts.res;

    Tm3(i) = mean(Ts);
    Tstd3(i) = std(Ts);

    sigs3(i) = (mean(Ts.^2)/(mean(Ts)^2)) - 1;

end


Tm1 = transpose(Tm1);
Tstd1 = transpose(Tstd1);


y1 = [Tm1(1:end-1) - Tstd1(1:end-1); Tm1(2:end) - Tstd1(2:end); flipud(Tm1(2:end) + Tstd1(2:end)); flipud(Tm1(1:end-1) + Tstd1(1:end-1))];

x1 = [dists(1:end-1); dists(2:end); dists(2:end); dists(1:end-1)];


Tm2 = transpose(Tm2);
Tstd2 = transpose(Tstd2);

y2 = [Tm2(1:end-1) - Tstd2(1:end-1); Tm2(2:end) - Tstd2(2:end); flipud(Tm2(2:end) + Tstd2(2:end)); flipud(Tm2(1:end-1) + Tstd2(1:end-1))];

x2 = [dists(1:end-1); dists(2:end); dists(2:end); dists(1:end-1)];


Tm3 = transpose(Tm3);
Tstd3 = transpose(Tstd3);

y3 = [Tm3(1:end-1) - Tstd3(1:end-1); Tm3(2:end) - Tstd3(2:end); flipud(Tm3(2:end) + Tstd3(2:end)); flipud(Tm3(1:end-1) + Tstd3(1:end-1))];

x3 = [dists(1:end-1); dists(2:end); dists(2:end); dists(1:end-1)];

hold on;
fill(x1, y1, 'red','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');
% errorbar(dists, Tm1, Tstd1,'LineWidth', 1.4, 'DisplayName', '$r_d=0.1m$');


% fill(x2, y2, 'red','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');
% plot(dists, Tm2, 'LineWidth', 1.4, 'DisplayName', '$r_d=0.15m$');
% errorbar(dists, Tm2, Tstd2,'LineWidth', 1.4, 'DisplayName', '$r_d=0.15m$');

fill(x3, y3, 'blue','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');



plot(dists, Tm3, 'v-','LineWidth', 1.4, 'DisplayName', '$r_d=0.2m$');
plot(dists, Tm1, 's-','LineWidth', 1.4, 'DisplayName', '$r_d=0.1m$');


% errorbar(dists, Tm3, Tstd3,'LineWidth', 1.4, 'DisplayName', '$r_d=0.2m$');

% errorbar(par(1:2:end), snr_m(1:2:end), snr_std(1:2:end),'-','LineWidth',1.2, 'CapSize',4, 'DisplayName','\Delta =10');



% errorbar(dists, Tm, Tstd, 'DisplayName', 'r_d=0.15m')





legend('Interpreter', 'latex')
xlabel('$L[m]$', 'Interpreter', 'latex');
ylabel('$T$', 'Interpreter', 'latex');
% savefigures('phasesc_T')
