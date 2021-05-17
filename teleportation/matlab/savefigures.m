function [] = savefigures(figure_name)

myfigure = gcf;
% pbaspect([1 0.618 1])

% myfigure.Position = [308,443,560,330];
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLineLineWidth', 1.2);

figure_name_full = ['fig_', figure_name];
saveas(myfigure, figure_name_full)
print(figure_name_full,'-dpdf','-r300',myfigure)

end