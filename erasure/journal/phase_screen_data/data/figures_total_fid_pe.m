clear all;
close all;
set(0,'defaultTextInterpreter','latex');

% newcolors = [0.9290, 0.6940, 0.1250
%              0.8500, 0.3250, 0.0980
%              0, 0.4470, 0.7410
%              0.25 0.80 0.54];
%          
% colororder(newcolors)


% PD
N = 10000;

rec = '0.2';
dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000];


% rec = '0.15';
% dists = [1000 1200 1400 1600 1800 2000 2200 2400 2600];

% rec = '0.1';
% dists = [1000 1200 1400 1600 1800 2000 2200 2400];

sigma_coh = 10;
F_class = (1 + (1/sigma_coh))/(2 + (1/sigma_coh));


data = load(['Fscatter_erasure_phasesc_rec=', rec,'_V3.mat']);
Fm_erasure3 = data.Fm;
Fmean_erasure3 = mean(transpose(Fm_erasure3));
Fstd_erasure3 = std(transpose(Fm_erasure3));

data = load(['Fscatter_erasure_phasesc_rec=', rec, '_V10.mat']);
Fm_erasure10 = data.Fm;
Fmean_erasure10 = mean(transpose(Fm_erasure10));
Fstd_erasure10 = std(transpose(Fm_erasure10));

data = load(['Fscatter_phasesc_rec=', rec,'_V3.mat']);
Fm3 = data.Fm;
Fmean3 = mean(transpose(Fm3));
Fstd3 = std(transpose(Fm3));

data = load(['Fscatter_phasesc_rec=', rec, '_V10.mat']);
Fm10 = data.Fm;
Fmean10 = mean(transpose(Fm10));
Fstd10 = std(transpose(Fm10));

data = load(['Fscatter_phasesc_rec=', rec, '_dir.mat']);
Fm_dir = data.Fm_dir;
Fmean_dir = mean(transpose(Fm_dir));
Fstd_dir = std(transpose(Fm_dir));
 
F0 = 0.;

pe1 = 0.1;
F_total1 = (1-pe1)^3*Fmean10 + (pe1*(1-pe1)^2)*(2*Fmean_erasure10 + Fmean10) + 3*(1-pe1)*pe1^2*F0 + pe1^3*F0;
F_total_dir1 = (1-pe1)*Fmean_dir + pe1*F0;

pe2= 0.2;
F_total2 = (1-pe2)^3*Fmean10 + (pe2*(1-pe2)^2)*(2*Fmean_erasure10 + Fmean10) + 3*(1-pe2)*pe2^2*F0 + pe2^3*F0;
F_total_dir2 = (1-pe2)*Fmean_dir + pe2*F0;

pe3 = 0.3;
F_total3 = (1-pe3)^3*Fmean10 + (pe3*(1-pe3)^2)*(2*Fmean_erasure10 + Fmean10) + 3*(1-pe3)*pe3^2*F0 + pe3^3*F0;
F_total_dir3 = (1-pe3)*Fmean_dir + pe3*F0;

pe4 = 0.4;
F_total4 = (1-pe4)^3*Fmean10 + (pe4*(1-pe4)^2)*(2*Fmean_erasure10 + Fmean10) + 3*(1-pe4)*pe4^2*F0 + pe4^3*F0;
F_total_dir4 = (1-pe4)*Fmean_dir + pe4*F0;


% data3 = load(['Fscatter_erasure_phasesc_rec=', rec, '_V6.mat']);
% Fm3 = data3.Fm;
% Fmean3 = mean(transpose(Fm3));
% Fstd3 = std(transpose(Fm3));


% data = load('F_pd_std.mat');
% Fstd = data.Fstd;



% y = [Fmean(1:end-1) - Fstd(1:end-1); Fmean(2:end) - Fstd(2:end); flipud(Fmean(2:end) + Fstd(2:end)); flipud(Fmean(1:end-1) + Fstd(1:end-1))];
% 
% y2 = [Fmean2(1:end-1) - Fstd2(1:end-1); Fmean2(2:end) - Fstd2(2:end); flipud(Fmean2(2:end) + Fstd2(2:end)); flipud(Fmean2(1:end-1) + Fstd2(1:end-1))];
% 
% x = [dists(1:end-1); dists(2:end); dists(2:end); dists(1:end-1)];
% 
% y_dir = [Fmean_dir(1:end-1) - Fstd_dir(1:end-1); Fmean_dir(2:end) - Fstd_dir(2:end); flipud(Fmean_dir(2:end) + Fstd_dir(2:end)); flipud(Fmean_dir(1:end-1) + Fstd_dir(1:end-1))];

% data = load('F_pd_std_dir.mat');
% Fstd_dir = data.Fstd_dir;

class = ones(size(dists))*F_class;

set(gca,'fontname','times') 
% subplot(2, 2, [1 3])
% figure;
hold on;

% yyaxis left
% title('Plots with Different y-Scales')
xlabel('$L[\mathrm{m}]$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');


% fill(x, y_dir, 'black','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');
% 
% fill(x, y, 'blue','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');
% 
% fill(x, y2, 'red','LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off');

color1 = [0, 0.4470, 0.7410];

plot(dists, F_total1, '-', 'Color', color1, 'LineWidth', 1.4, 'DisplayName', '$p_\mathrm{e}=0.1$');

plot(dists, F_total_dir1, '--', 'Color', color1, 'LineWidth', 1.4, 'HandleVisibility','off');

color2 = [0.8500, 0.3250, 0.0980];

plot(dists, F_total2, '-', 'Color', color2, 'LineWidth', 1.4, 'DisplayName', '$p_\mathrm{e}=0.2$');

plot(dists, F_total_dir2, '--', 'Color', color2, 'LineWidth', 1.4, 'HandleVisibility','off');

color3 = [0.9290, 0.6940, 0.1250];

plot(dists, F_total3, '-', 'Color', color3, 'LineWidth', 1.4, 'DisplayName', '$p_\mathrm{e}=0.3$');

plot(dists, F_total_dir3, '--', 'Color', color3, 'LineWidth', 1.4, 'HandleVisibility','off');

color4 = [0.4940 0.1840 0.5560];

plot(dists, F_total4, '-', 'Color', color4, 'LineWidth', 1.4, 'DisplayName', '$p_\mathrm{e}=0.4$');

plot(dists, F_total_dir4, '--', 'Color', color4, 'LineWidth', 1.4, 'HandleVisibility','off');

% xlim([par(1) par(end)])
ylim([0.5 0.95]);
% xlim([1000 1910]);
% legend('fontname','times', 'Location','northeast','Interpreter', 'latex');
% text(1020, .98, ['$r_\mathrm{d}=', rec, 'm$']);

txt = ['$r_\mathrm{d} =', rec, '$ m'];
% text(1750, .80, txt);
text(2700, .93, txt);

savefigures('total2');



% yyaxis right
% ylabel('$\%> F_\mathrm{class}$', 'Interpreter', 'latex');
% 
% plot(par, Ncl, 'LineWidth', 1.7, 'DisplayName', 'Protocol');
% plot(par, Ncl_dir, 'LineWidth', 1.7, 'DisplayName', 'Direct');


