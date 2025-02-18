clear;
set(0,'defaultTextInterpreter','latex');
% Fidelity parameters
V = 10;
coh_sigma = 10;
epsilon = 0.05;
eta = sqrt(.9);

% Distribution parameters
% Discrete error rate
pd = 0.;
% Discrete error transmissivity
T0 = 0.;
% Gaussian mean 
mu = 0.8;
% Gaussian sigma
sigma = 0.2; 
% figure;
set(gca,'fontname','times') 
hold on;
xlabel('$g$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');

% Sample size
N = 5;
Ts1 = distribution(pd, T0, mu, sigma, N);
Ts2 = distribution(pd, T0, mu, sigma, N);
Ts3 = distribution(pd, T0, mu, sigma, N);

Ts1 = [.87 0. 0.75 0.88 0.66];
Ts2 = [.92 0.89 0.72 0. 0.42];
Ts3 = [0. 0.94 0.76 0.83 0.71];
% 
% Ts1 = [0.8 0.6 0.7 .3];
% Ts2 = [0.8 0.2 0.7 0.8];
% Ts3 = [0.4 0.7 0.7 0.6];

N = length(Ts1);
Fmax = zeros(N,1);
gmax = zeros(N, 1);
gs = -1.12:.01:1.12;

for j = 1:N
    T1 = Ts1(j);
    T2 = Ts2(j);
    T3 = Ts3(j);

    disp([num2str((fid_tmsv_dir_eq(T1, epsilon, coh_sigma) ...
            +fid_tmsv_dir_eq(T2, epsilon, coh_sigma) ...
            +fid_tmsv_dir_eq(T3, epsilon, coh_sigma))/3)]);

    Fs = fid_tmsv_gen_loss_eq(V, gs, coh_sigma, T1, T2, T3, epsilon, eta);
    


% %     plot(gs, Fs);
    plot(gs, Fs, 'LineWidth', 1.)
    [Fmax(j), gmax(j)]= max(Fs);
    Ttext = ['[' num2str(round(T1,2)) ' , ' num2str(round(T2,2))  ' , ' num2str(round(T3,2)) ']' ];
    text(gs(gmax(j)) - 0.14, Fmax(j) + 0.03, Ttext);
end

scatter(gs(gmax), Fmax, 'filled', 'MarkerFaceColor', 'red');
xlim([-1.11, 1.11])
ylim([0.3, 1])
savefigures('gs')


