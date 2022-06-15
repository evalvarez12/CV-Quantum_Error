
V = 20;
sigma = 10;
epsilon = 0;

sigmaT = 0.6;
meanT = 0.5;

N = 2e6;
Ts = normrnd(meanT, sigmaT, N, 3);
Ts = max(Ts, 0);
Ts(Ts>1) = 1;

Fmax = zeros(N,1);
gmax = zeros(N, 1);

gs = -1:.01:1;



for i = 1:N
    T1 = Ts(i, 1);
    T2 = Ts(i, 2);
    T3 = Ts(i, 3);

    Fs = fid_tmsv_gen_loss_eq(V, gs, sigma, T1, T2, T3, epsilon);
    [Fmax(i), gmax(i)]= max(Fs);


end


mean(Fmax)

