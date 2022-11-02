clear;

V = 10;
sigma = 10;

g = -1;

% T = 0.9636 - d = 0.4032
% T = 1 - d = 0.4049
% T = 0.6 - d = 0.3941
% T = 0.634 - d =  0.3944
% T = 0.6727 - d = 0.3949

T1 = 0;
T2 = 1;
T3 = 1;

da = 30;
di = 0.01;

d = 0.4032;


xbar_c = da;
pbar_c = di;

bar_d = 0.05;
bar_N = 2;


xbar_ini = xbar_c - (bar_N/2)*bar_d;
xbar_end = xbar_c + (bar_N/2)*bar_d;
pbar_ini = pbar_c - (bar_N/2)*bar_d;
pbar_end = pbar_c + (bar_N/2)*bar_d;

xbars = xbar_ini:bar_d:xbar_end;
pbars = pbar_ini:bar_d:pbar_end;

ps = zeros(bar_N, bar_N);


for i=1:bar_N
    xbar = xbars(i);

    for j=1:bar_N
        i
        j
        pbar = pbars(j);
    
        ps(i, j) = prob_dist_sb(xbar, pbar, da, di, T1, T2, T3, V, d, sigma);
    end
end






