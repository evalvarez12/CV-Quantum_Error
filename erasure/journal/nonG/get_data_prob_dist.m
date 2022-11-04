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

d = 0.4032;


da = 30;
di = 0.01;

delta = 5;
N = 5;

xbars = linspace(da-delta, da+delta,N);
pbars = linspace(di-delta, di+delta,N);

ps = zeros(N^2, 1);

[x, y] = meshgrid(xbars, pbars);
xf = x(:);
yf = y(:);



parfor i=1:N^2
    i
    xbar = xf(i);
    pbar = yf(i);
    
    ps(i) = prob_dist_sb(xbar, pbar, da, di, T1, T2, T3, V, d, sigma);
   
end


ps = reshape(ps, N, N);





