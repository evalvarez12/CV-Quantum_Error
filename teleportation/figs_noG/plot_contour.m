clear all;
close all;


data_py = load('data-py-sig10');
data_as = load('data-as-sig10');

tmsv = zeros(30, 30);
sb = zeros(30, 30);
for i = 1:30
    for j = 1:30
        tmsv(i,j) = data_py.data(1,i,j);
        sb(i,j) = data_py.data(2,i,j);
    end
end


as = data_as.results;


ts = linspace(0.0001, 1, 30);
es = logspace(-2.2, -0.9, 30);



mask = as > sb;

top = sb;
top(mask) = as(mask);

color_mask = zeros(30, 30);
color_mask(mask) = 1;
color_mask(~mask) =  0;
color_mask(as < 0.5) = -1;



figure 
contourf(es, ts, top);



figure; hAxes = gca;

hold on
set(gca,'FontSize',10)

[C,h] = contour(es, ts, top, 'LevelList', [0.6, 0.65, 0.7, 0.75, 0.8, 0.85], 'EdgeColor', 'k');   
% clabel(C,h)
contP = get(h,'Parent');

% set(findobj(gca,'Type','patch','UserData',2),'EdgeColor',[0 0 0])

X = contP.XLim;
Y = contP.YLim;

% ax1.Visible = 'off';
colormap('hot')

% mnew = imresize(color_mask,30,'bilinear');


img = imagesc(color_mask,'XData',X,'YData',Y);


colormap( hAxes , [1 0 0; 0 0 1 ;0 1 0] );
img.AlphaData = 0.2;


% [C,h] = contour(ts, es, top);   
clabel(C,h,'FontSize', 10);


set(gca, 'XScale', 'log');
xlim([es(1), es(end)]);

ylabel('$T$', 'Interpreter', 'latex');
xlabel('$\epsilon$', 'Interpreter', 'latex');



% text markers 
txt_sat1 = ('s1');
txt_sat2 = ('s2');
txt_sat3 = ('s3');
txt_sat4 = ('s4');

txt_fib11 = ('f11');
txt_fib12 = ('f12');
txt_fib13 = ('f13');
txt_fib14 = ('f14');

txt_fib21 = ('f21');
txt_fib22 = ('f22');
txt_fib23 = ('f23');
txt_fib24 = ('f24');

text(0.008, .158, txt_fib11);
text(0.008, .0631, txt_fib12);
text(0.026, .158, txt_fib13);
text(0.026, .0631, txt_fib14);

text(0.01, .0251, txt_fib21);
text(0.01, .004, txt_fib22);
text(0.028, .0251, txt_fib23);
text(0.028, .004, txt_fib24);

text(0.0135, .063, txt_sat1);
text(0.0135, .01, txt_sat2);
text(0.0144, .063, txt_sat3);
text(0.0144, .01, txt_sat4);

