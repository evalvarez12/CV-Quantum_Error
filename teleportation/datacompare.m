clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Fixed loss
mydata = load('matlab/tmsv.mat');
mydata = mydata.data;
mydata = transpose(mydata);

mgdata = load('MGsCFcode/results.mat');
mgdata = mgdata.results;
mgdata = transpose(mgdata);

mgtmsv = mgdata(1,:);
mgps = mgdata(2,:);
mgas = mgdata(3,:);


fig = figure;

hold on

zs = linspace(0, 1, 15);

plot(zs, mydata, 'o', 'DisplayName', 'My TMSV', 'color', [0, 0, 1]);
plot(zs, mgtmsv, '-', 'DisplayName', 'Mingjian TMSV', 'color', 	[1, 0, 0]);
plot(zs, mgas, '-', 'DisplayName', 'Mingjian PA-PS', 'color', [0, 0.5, 0]);
% plot(zs, tele4(1:1:end), '-', 'DisplayName', 'Teleportation $\sigma = 25$', 'color', [0.4940, 0.1840, 0.5560] );
% 
% ylim([0.46 1]);
% xlim([0 13]);

ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex')
xlabel('Transmissivity', 'Interpreter', 'latex')
legend('Location','northwest');

% grid on;
 set(gca,'Box','on');

savefigures('fig_compare');