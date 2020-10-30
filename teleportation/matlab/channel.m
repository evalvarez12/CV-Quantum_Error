clear all;
close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Histograms
t_ext = 10^(-1/10) * exp(-0.7);

data = load('../../../Laser_propagation/data/TELE_DOWN_I_r=1_z=0_1024_10000');
data = data.res .* t_ext;

figure
histogram(data, 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
xline(mean(data), 'k--', 'LineWidth', 1.25);
xlim([min(data) max(data)]);

ylabel('PDF', 'Interpreter', 'latex');
xlabel('$T$', 'Interpreter', 'latex');

txt1 = ('Downlink');
txt2 = ('$L_0=\infty$ $l_0=0$');
text(0.205 , 650, txt1);
text(0.205 , 550, txt2);

savefigures('hist1');


data = load('../../../Laser_propagation/data/TELE_UP_I_r=1_z=0_1024_10000');
data = data.res .* t_ext;

figure
histogram(data, 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6);
xline(mean(data), 'k--', 'LineWidth', 1.25);
xlim([min(data) max(data)]);

ylabel('PDF', 'Interpreter', 'latex');
xlabel('$T$', 'Interpreter', 'latex');

txt1 = ('Uplink');
txt2 = ('$L_0=$Coulman Vernin $l_0=\delta L_0$');
text(0.06 , 21, txt1);
text(0.06 , 17, txt2);

savefigures('hist2');


data = load('../../../Laser_propagation/data/TELE__L0inf_UP_I_r=1_z=0_1024_10000');
data = data.res .* t_ext;

figure
histogram(data, 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.6, 'NumBins', 35);
xline(mean(data), 'k--', 'LineWidth', 1.25);
xlim([min(data) max(data)]);

ylabel('PDF', 'Interpreter', 'latex');
xlabel('$T$', 'Interpreter', 'latex');

txt1 = ('Uplink');
txt2 = ('$L_0=\infty$ $l_0=0$');
text(0.12 , 12, txt1);
text(0.12 , 10, txt2);
savefigures('hist3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Channel properties
data = load('channel.mat');
data = data.data;

zs = data(1,:);
d_t = data(2,:);
d_noise = data(3,:);
u_t = data(4,:);
u_noise = data(5,:);

figure
% title('Satellite communications channels properties')
hold on

yyaxis left
plot(zs, d_t, 'o-', 'DisplayName', 'Teleportation');
plot(zs, u_t, 'o--', 'DisplayName', 'Direct');
ylabel('$T_f$ [dB]', 'Interpreter', 'latex')


yyaxis right
plot(zs, d_noise, 'o-', 'DisplayName', 'Teleportation');
plot(zs, u_noise, 'o--', 'DisplayName', 'Direct');
ylabel('$\epsilon_f$', 'Interpreter', 'latex')
set ( gca, 'ydir', 'reverse' )
ylim([-0.005 0.085]);

txt1 = ('$1.1 \times 10^{-5}$');
text(5 , 0.01, txt1);

txt2 = ('$1.2 \times 10^{-4}$');
text(30 , 0.01, txt2);

legend
xlabel('$\zeta$ [deg]', 'Interpreter', 'latex')
% grid on;
 set(gca,'Box','on');

savefigures('channel');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Fixed loss
data = load('fixed.mat');
data = data.data;

zs = data(1,:);
tele1 = data(2,:);
tele2 = data(3,:);
tele3 = data(4,:);
tele4 = data(5,:);
dir1 = data(6,:);
dir2 = data(7,:);
dir3 = data(8,:);
dir4 = data(9,:);

classical = ones(1, length(tele1)) * .5;

fig = figure;
% title('Fixed channel fidelities')

hold on

plot(zs(1:1:end), tele1(1:1:end), '-', 'DisplayName', 'Teleportation $\sigma = 2$', 'color', [0, 0, 1]);
plot(zs(1:1:end), tele2(1:1:end), '-', 'DisplayName', 'Teleportation $\sigma = 4$', 'color', 	[1, 0, 0]);
plot(zs(1:1:end), tele3(1:1:end), '-', 'DisplayName', 'Teleportation $\sigma = 10$', 'color', [0, 0.5, 0]);
plot(zs(1:1:end), tele4(1:1:end), '-', 'DisplayName', 'Teleportation $\sigma = 25$', 'color', [0.4940, 0.1840, 0.5560] );
plot(zs(1:1:end), dir1(1:1:end), '--', 'DisplayName', 'Direct $\sigma = 2$', 'color', [0, 0, 1]);
plot(zs(1:1:end), dir2(1:1:end), '--', 'DisplayName', 'Direct $\sigma = 4$', 'color', 	[1, 0, 0]);
plot(zs(1:1:end), dir3(1:1:end), '--', 'DisplayName', 'Direct $\sigma = 10$', 'color', [0, 0.5, 0]);
plot(zs(1:1:end), dir4(1:1:end), '--', 'DisplayName', 'Direct $\sigma = 25$', 'color', [0.4940, 0.1840, 0.5560]);

plot(zs, classical, 'r--', 'HandleVisibility','off');

ylim([0.46 1]);
xlim([0 13]);

ylabel('$\bar{F}$', 'Interpreter', 'latex')
xlabel('Fixed loss [dB]', 'Interpreter', 'latex')
legend;

txt1 = ('classical limit');
text(0.1 , 0.52, txt1, 'Color', 'r');

% grid on;
 set(gca,'Box','on');

savefigures('fixed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Fidelities
data = load('fidelities.mat');
data = data.data;

zs = data(1,:);
tele1 = data(2,:);
tele2 = data(3,:);
tele3 = data(4,:);
tele4 = data(5,:);
dir1 = data(6,:);
dir2 = data(7,:);
dir3 = data(8,:);

classical = ones(1, length(tele1)) * .5;

fig = figure;
% title('Satellite-to-ground fidelities')

hold on

plot(zs, tele1, 'o-', 'DisplayName', 'Teleportation $\sigma = 2$', 'color', [0, 0, 1]);
plot(zs, tele2, 'o-', 'DisplayName', 'Teleportation $\sigma = 4$', 'color', 	[1, 0, 0]);
plot(zs, tele3, 'o-', 'DisplayName', 'Teleportation $\sigma = 10$', 'color', [0, 0.5, 0]);
plot(zs, tele4, 'o-', 'DisplayName', 'Teleportation $\sigma = 25$', 'color', [0.4940, 0.1840, 0.5560] );
plot(zs, dir1, 'v--', 'DisplayName', 'Direct $\sigma = 2$', 'color', [0, 0, 1]);
plot(zs, dir2, 'v--', 'DisplayName', 'Direct $\sigma = 4$', 'color', 	[1, 0, 0]);
% plot(zs, dir3, 'v--', 'DisplayName', 'Direct $\sigma = 10$', 'color', 	[0, 0.5, 0]);

plot(zs, classical, 'r--', 'HandleVisibility','off');

ylim([0.35 0.72]);

ylabel('$\bar{F}$', 'Interpreter', 'latex')
xlabel('$\zeta$ [deg]', 'Interpreter', 'latex')
legend('Location', 'northeast');

txt1 = ('classical limit');
text(0.1 , 0.48, txt1, 'Color', 'r');

% grid on;
 set(gca,'Box','on');

savefigures('fidelities');

