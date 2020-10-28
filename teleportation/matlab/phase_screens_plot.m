h0 = 0;
h1 = 20e3;
H = 500e3;
k = 2 * pi / 1.55e-6;

i = 0:.1:h1;

f = @(x) HVmodel(x);
hs = solve_boundaries(f, h0, h1, 18);
hs = hs(1:end-1);

hs3 = (0:20)*h1/18;

hs2 = solve_equal_rytov(k, h0, h1, 0.2);
disp(['num ph: ', num2str(length(hs2))])

figure('DefaultAxesFontSize',14)

hold on

plot(HVmodel(i), i);
% plot(log10(HVmodel(hs)), hs, 'ro');
plot(HVmodel(hs3), hs3, 'r*');
plot(HVmodel(hs2), hs2, 'bo');

leg = legend('H-V${}_{5/6}$ Model', 'Uniform', 'Equal Rytov');
set(leg, 'Interpreter', 'latex');
% set(leg1,'FontSize',17);

ylim([-100 h1]);
xlim([7e-19 1.8e-14]);

set(gca, 'XScale', 'log');
set(gca,'fontname','times')  % Set it to times

xlabel('$C_n^2(h)$', 'Interpreter', 'latex')
ylabel('$h$', 'Interpreter', 'latex')

grid on

savefigures('phase_screens');