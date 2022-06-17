N = 1e5;
pd = 0.5;
T0 = 0.;

sigma = 0;
mean = 1;

s = rand(N,1);
s(s>pd) = 1;
s(s ~= 1) = T0;


ts = normrnd(mean, sigma, N, 1);
ts = ts .* s;

ts = max(ts, 0);


extra = length(ts(ts>1));
i = randi([1 N-extra], extra, 1);
support = ts(ts <= 1);

ts(ts > 1) = support(i);


% histogram(s);
histogram(ts, 100);


% X = 0:.01:1;
% y1 = betapdf(X,10,9);
% y2 = betapdf(X,10,8);
% y3 = betapdf(X,10,4);
% y4 = betapdf(X,10,2);
% y5 = betapdf(X,10,1);
% 
% figure
% plot(X,y1,'Color','r','LineWidth',2)
% hold on
% plot(X,y2,'LineStyle','-.','Color','b','LineWidth',2)
% plot(X,y3,'LineStyle',':','Color','g','LineWidth',2)
% plot(X,y4,'LineStyle','--','Color','k','LineWidth',2)
% plot(X,y5,'LineStyle','--','Color','y','LineWidth',2)
% 
% legend({'a = b = 0.75','a = b = 1','a = b = 4'},'Location','NorthEast');
% hold off

% 
% X = 0:.01:1;
% y1 = lognpdf(X,0,1);
% y2 = lognpdf(X,4,1);
% y3 = lognpdf(X,1,3);
% y4 = lognpdf(X,1,2);
% y5 = lognpdf(X,1,1);;
% 
% figure
% plot(X,y1,'Color','r','LineWidth',2)
% hold on
% plot(X,y2,'LineStyle','-.','Color','b','LineWidth',2)
% plot(X,y3,'LineStyle',':','Color','g','LineWidth',2)
% plot(X,y4,'LineStyle','--','Color','k','LineWidth',2)
% plot(X,y5,'LineStyle','--','Color','y','LineWidth',2)
% 
% legend({'a = b = 0.75','a = b = 1','a = b = 4'},'Location','NorthEast');
% hold off


% clear
% x = 0: 0.001 : 1;
% x = x';
% y1 = gampdf(x,1.5,2);
% y1 = y1 ./ sum(y1);
% maxx1 = x(y1==max(y1)),expx1 = 1.5/2
% y2 = gampdf(x,2.0,2);
% y2 = y2 ./ sum(y2);
% maxx2 = x(y2==max(y2)),expx2 = 2/2
% y3 = gampdf(x,1.5+2,2);
% y3 = y3 ./ sum(y3);
% maxx3 = x(y3==max(y3)),expx3 = (1.5+2)/2
% 
% xc = 2*x;
% yc = conv(y1,y2)*(x(2)-x(1));
% yc = interp1(xc,yc(1:2:end),x,...
% 'linear','extrap'); xc = x;
% yc = yc./sum(yc);
% maxxc = x(yc==max(yc))
% 
% figure2 = figure('Color',[1 1 1],...
%    'Position',[50 200 600 300]);
% patch(x,y1,[0.1 0.3 0.8]), hold on
% % patch(x,y2,[0.1 0.8 0.1])
% % patch(x,yc,[1 0 0])
% alpha(0.3)
% title('Convolution of Distributions')
% grid, hold off