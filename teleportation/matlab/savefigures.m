function [] = savefigures(figure_name)

myfigure = gcf;
myfigure.Position = [308,443,560,330];
figure_name_full = ['fig_', figure_name];
saveas(myfigure, figure_name_full)
print(figure_name_full,'-depsc','-r300',myfigure)

end