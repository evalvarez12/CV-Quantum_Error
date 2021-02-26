clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Channel properties
data = load('channel_assym.mat');
data = data.data;

zs = data(1,:);
d_t = data(2,:);
d_noise = data(3,:);
u_t = data(4,:);
u_noise = sort(data(5,:));

figure
% title('Satellite communications channels properties')
hold on

yyaxis left
plot(zs, d_t, 'o-', 'DisplayName', 'Teleportation');
plot(zs, u_t, 'o--', 'DisplayName', 'Direct');
ylabel('$T_f$ [dB]', 'Interpreter', 'latex')
ylim([11 27]);


yyaxis right
plot(zs, d_noise, 'o-', 'DisplayName', 'Teleportation');
plot(zs, u_noise, 'o--', 'DisplayName', 'Direct');
ylabel('$\epsilon_f$', 'Interpreter', 'latex')
set ( gca, 'ydir', 'reverse' )
ylim([0 0.205]);

txt1 = ('$7.1 \times 10^{-3}$');
% text(48 , 0.015, txt1);

txt2 = ('$7.0 \times 10^{-3}$');
% text(2 , 0.015, txt2);

legend('Position', [.16 .51 .2 .2]);
xlabel('$\zeta$ [deg]', 'Interpreter', 'latex')
% grid on;
 set(gca,'Box','on');

savefigures('channel_assym');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Fidelities
data = load('fidelities_assym.mat');
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


fill([0,60,60,0], [0,0,0.5,0.5],  [1,.75,.75] ,'Edgecolor', 'none', 'handlevisibility', 'off')


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

ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex')
xlabel('$\zeta$ [deg]', 'Interpreter', 'latex')
legend('Location', 'northeast');

txt1 = ('classical limit');
text(0.8 , 0.48, txt1, 'Color', 'r');

% grid on;
 set(gca,'Box','on');
legend('NumColumns', 2);
 
savefigures('fidelities_assym');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Fidelities optimized direct
data = load('fidelities_assym_opt_direct.mat');
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


fill([0,60,60,0], [0,0,0.5,0.5],  [1,.75,.75] ,'Edgecolor', 'none', 'handlevisibility', 'off')


hold on

plot(zs, tele1, 'o-', 'DisplayName', 'Teleportation $\sigma = 2$', 'color', [0, 0, 1]);
plot(zs, tele2, 'o-', 'DisplayName', 'Teleportation $\sigma = 4$', 'color', 	[1, 0, 0]);
plot(zs, tele3, 'o-', 'DisplayName', 'Teleportation $\sigma = 10$', 'color', [0, 0.5, 0]);
plot(zs, tele4, 'o-', 'DisplayName', 'Teleportation $\sigma = 25$', 'color', [0.4940, 0.1840, 0.5560] );
plot(zs, dir1, 'v--', 'DisplayName', 'Direct $\sigma = 2$', 'color', [0, 0, 1]);
plot(zs, dir2, 'v--', 'DisplayName', 'Direct $\sigma = 4$', 'color', 	[1, 0, 0]);
plot(zs, dir3, 'v--', 'DisplayName', 'Direct $\sigma = 10$', 'color', 	[0, 0.5, 0]);

plot(zs, classical, 'r--', 'HandleVisibility','off');

ylim([0.07 0.65]);

ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex')
xlabel('$\zeta$ [deg]', 'Interpreter', 'latex')
legend('Location', 'southeast');

txt1 = ('classical limit');
% text(0.8 , 0.48, txt1, 'Color', 'r');

% grid on;
 set(gca,'Box','on');
legend('NumColumns', 2);
 
savefigures('fidelities_assym_opt_direct');