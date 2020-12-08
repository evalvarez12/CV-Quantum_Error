clc;clear

OptOption = optimoptions(@fmincon, 'FunctionTolerance', 1e-30,'StepTolerance', 1e-20, 'Display','off');

% 'algorithm','sqp'
e = 1e-4;

rdB = [5, 5, 10];
sig = [1, sqrt(10), sqrt(10)];

rnum = 3;

r = rdB/8.67;

tnum = 32;
tdB = linspace(1, 15, tnum);

F_ps = zeros(rnum, tnum);
F_pa = F_ps;
F_pc = F_ps;
F_as = F_ps;
F_sa = F_ps;
Fepr = F_ps;

opt_ps = F_ps;
opt_pa = F_ps;
opt_pc = F_ps;
opt_sa = F_ps;
opt_as = F_ps;


gmax = 10;
gmin = 0.1;
gini = 0.1;

% gmax = 1;
% gmin = 1;
% gini = 1;

Tmax = 0.999;
Tmin = 0.001;
Tini = 0.99;

tic
for i = 1:rnum
    disp(i);
    parfor j = 1:tnum
        l = tanh(r(i));
        e_eff = e * 2*cosh(r(i));
        t = 10.^(-tdB(j)/10);
        
        fun_ep = @(par) -coh_loss(par(1), l, t, e_eff, 'epr', par(2), sig(i));
        fun_ps = @(par) -coh_loss(par(1), l, t, e_eff, 'ps', par(2),  sig(i));
        fun_pa = @(par) -coh_loss(par(1), l, t, e_eff, 'pa', par(2),  sig(i));
        fun_pc = @(par) -coh_loss(par(1), l, t, e_eff, 'pc', par(2),  sig(i));
        fun_as = @(par) -coh_loss(par(1), l, t, e_eff, 'as', par(2),  sig(i));
        fun_sa = @(par) -coh_loss(par(1), l, t, e_eff, 'sa', par(2),  sig(i));
        
        [~, Fepr(i,j)] = fmincon(fun_ep, [Tini, gini], [],[],[],[], [Tmin, gmin], [Tmax, gmax], [], OptOption);
        [~, F_ps(i,j)] = fmincon(fun_ps, [Tini, gini], [],[],[],[], [Tmin, gmin], [Tmax, gmax], [], OptOption);
        [~, F_pa(i,j)] = fmincon(fun_pa, [Tini, gini], [],[],[],[], [Tmin, gmin], [Tmax, gmax], [], OptOption);
        [~, F_pc(i,j)] = fmincon(fun_pc, [Tini, gini], [],[],[],[], [Tmin, gmin], [Tmax, gmax], [], OptOption);
        [~, F_as(i,j)] = fmincon(fun_as, [Tini, gini], [],[],[],[], [Tmin, gmin], [Tmax, gmax], [], OptOption);
        [~, F_sa(i,j)] = fmincon(fun_sa, [Tini, gini], [],[],[],[], [Tmin, gmin], [Tmax, gmax], [], OptOption);
        
    end
end
toc

Fepr = -Fepr;
F_ps = -F_ps;
F_pa = -F_pa;
F_pc = -F_pc;
F_as = -F_as;
F_sa = -F_sa;

%%
% save('coh_fmax_s10.mat')
FontSize = 12;

PRINT = 1;
lw = 1.25;

MarkerIndices = 1:4:tnum;
%
rindex = 1;
figure
fill([0,15,15,0], [0,0,0.5,0.5],  [1,.75,.75] ,'Edgecolor', 'none', 'handlevisibility', 'off')
hold on
plot(tdB, F_ps(rindex,:), 'bo-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_pa(rindex,:), 'ro-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_pc(rindex,:), 'go-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_sa(rindex,:), 'b*-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_as(rindex,:), 'r*-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, Fepr(rindex,:),'kv-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, .5*ones(1, 32), 'r--', 'linewidth', lw, 'handlevisibility', 'off')


legend('PS','PA','PC','PS-PA', 'PA-PS', 'TMSV', 'numcolumns', 2, ...
    'location', 'northwest', 'FontSize', FontSize, 'Position', ...
    [0.57797619047619,0.705252527853455,0.306349440687698,0.18969696709604])

ylabel('$\bar\mathcal{F}$', 'FontSize', FontSize)
xlabel('$T_{\mathrm{f}}\,\textrm{[dB]}$', 'FontSize', FontSize)

xlim([1,15])
ylim([0.1, 1])

myfigure = gcf;
myfigure.Position = [308,443,560,330];
pbaspect([1 .618 1])

% Create text
rtext = text(1,1,'classical limit');
rtext.Position = [1.33870967741935 0.164230769230769 0];
rtext.Color = [1 0 0];
figure1 = gcf;
% Create arrow
annotation(figure1,'arrow',[0.155357142857143 0.214285714285714],...
    [0.462121212121212 0.239393939393939],'Color',[1 0 0]);

rtext = text(1,1,'$r=5\, \rm{[dB]}, \sigma=1$');
rtext.Units = 'normalized';
rtext.Position = [0.258064516129032,0.912932073681128,0];

if PRINT
    figure_name = 'fig_cohr5dBs1';
    saveas(myfigure,figure_name)
    print(figure_name,'-depsc','-r300',myfigure)
end
%
rindex = 2;
figure
fill([0,15,15,0], [0,0,0.5,0.5],  [1,.75,.75] ,'Edgecolor', 'none', 'handlevisibility', 'off')
hold on
plot(tdB, F_ps(rindex,:), 'bo-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_pa(rindex,:), 'ro-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_pc(rindex,:), 'go-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_sa(rindex,:), 'b*-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_as(rindex,:), 'r*-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, Fepr(rindex,:),'kv-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, .5*ones(1, 32), 'r--', 'linewidth', lw, 'handlevisibility', 'off')

legend('PS','PA','PC','PS-PA', 'PA-PS', 'TMSV', 'numcolumns', 2, ...
    'location', 'northwest', 'FontSize', FontSize, 'Position', ...
    [0.57797619047619,0.705252527853455,0.306349440687698,0.18969696709604])

ylabel('$\bar\mathcal{F}$', 'FontSize', FontSize)
xlabel('$T_{\mathrm{f}}\,\textrm{[dB]}$', 'FontSize', FontSize)

xlim([1,15])
ylim([0.1, 1])

myfigure = gcf;
myfigure.Position = [308,443,560,330];
pbaspect([1 0.618 1])

rtext = text(1,1,'$r=5\,\rm{[dB]}, \sigma=10$');
rtext.Units = 'normalized';
rtext.Position = [0.258064516129032,0.912932073681128,0];

% Create text
rtext = text(1,1,'classical limit');
rtext.Position = [1.33870967741935 0.164230769230769 0];
rtext.Color = [1 0 0];
figure1 = gcf;
% Create arrow
annotation(figure1,'arrow',[0.155357142857143 0.214285714285714],...
    [0.462121212121212 0.239393939393939],'Color',[1 0 0]);

if PRINT
    figure_name = 'fig_cohr5dBs10';
    saveas(myfigure,figure_name)
    print(figure_name,'-depsc','-r300',myfigure)
end

%
rindex = 3;
figure
fill([0,15,15,0], [0,0,0.5,0.5],  [1,.75,.75] ,'Edgecolor', 'none', 'handlevisibility', 'off')
hold on
plot(tdB, F_ps(rindex,:), 'bo-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_pa(rindex,:), 'ro-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_pc(rindex,:), 'go-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_sa(rindex,:), 'b*-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, F_as(rindex,:), 'r*-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, Fepr(rindex,:),'kv-', 'linewidth', lw, 'MarkerIndices', MarkerIndices)
plot(tdB, .5*ones(1, 32), 'r--', 'linewidth', lw, 'handlevisibility', 'off')

legend('PS','PA','PC','PS-PA', 'PA-PS', 'TMSV', 'numcolumns', 2, ...
    'location', 'northwest', 'FontSize', FontSize, 'Position', ...
    [0.57797619047619,0.705252527853455,0.306349440687698,0.18969696709604])

ylabel('$\bar\mathcal{F}$', 'FontSize', FontSize)
xlabel('$T_{\mathrm{f}}\,\textrm{[dB]}$', 'FontSize', FontSize)

xlim([1,15])
ylim([0.1, 1])

myfigure = gcf;
myfigure.Position = [308,443,560,330];
pbaspect([1.618 1 1])

rtext = text(1,1,'$r=10\, \rm{[dB]}, \sigma=10$');
rtext.Units = 'normalized';
rtext.Position = [0.258064516129032,0.912932073681128,0];

% Create text
rtext = text(1,1,'classical limit');
rtext.Position = [1.33870967741935 0.164230769230769 0];
rtext.Color = [1 0 0];
figure1 = gcf;
% Create arrow
annotation(figure1,'arrow',[0.155357142857143 0.214285714285714],...
    [0.462121212121212 0.239393939393939],'Color',[1 0 0]);

if PRINT
    figure_name = 'fig_cohr10dBs10';
    saveas(myfigure,figure_name)
    print(figure_name,'-depsc','-r300',myfigure)
end


%%
% figure
% Fd = F_pas;
% surf(rdB, tdB, Fepr','LineStyle','--',...
%     'FaceAlpha',0.5,...
%     'FaceColor',[0.149019607843137 0.149019607843137 0.149019607843137]);
% hold on
% surf(rdB, tdB, Fd','LineStyle','none',...
%     'FaceColor',[0.074509803921569 0.623529411764706 1]);
% ylabel('$T$[dB]')
% xlabel('$r$[dB]')
% view([-129.317647058823 23.1116682678502]);
% grid('on');
% hold('off');