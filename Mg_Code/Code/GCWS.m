clc
clear
close all
warning('off')
sum_limit_n = 30;
sum_limit_m = 20;
% Calculate z(n,k,m,l)
z = zeros(sum_limit_n + 2, sum_limit_n + 2, sum_limit_m + 1, sum_limit_m + 1);
for n = 0:sum_limit_n + 2
    for k = 0:n
        for m = 0:sum_limit_m
            for l = 0:m
                z(n + 1, k + 1, m + 1, l + 1) = sqrt(nchoosek(n - k + l, l))...
                    * sqrt(nchoosek(k + m - l, k));
            end
        end
    end
end

th = 1e-20;
Tbs = 0.9;
Nbin = [1 2];
Tbin = 10.^(-(linspace(.1,150,100)/50));
B = .001;
A = 1.3;

ParaPho.B = B;
ParaPho.th = th;
ParaPho.Tbs = Tbs;
ParaPho.sum_limit_n = sum_limit_n;
ParaPho.sum_limit_m = sum_limit_m;
ParaPho.mode = 'GetMinT';


%% Fixed
tic
ParaPho.operation = 'No';
[RateNo, DistanceNo,~] = Operation(1, 1, ParaPho, A, z, 0);

ParaPho.operation = 'TxS';
[RateT, DistanceT,~] = Operation(1, 1, ParaPho, A, z, 0);

ParaPho.operation = 'RxS';
[RateR, DistanceR,~] = Operation(1, 1, ParaPho, A, z, 0);
toc

Attenuation = 0;%8.3746;
%semilogy(-50*log10(Tbin), -log2(1-Tbin),'k--','linewidth',1.5);hold on
figure(1)
semilogy(10.^((DistanceNo - Attenuation)/-50), RateNo ,'k--','linewidth',1.5);
hold on
semilogy(10.^((DistanceT - Attenuation)/-50), RateT,'b','linewidth',1.5);
semilogy(10.^((DistanceR - Attenuation)/-50), RateR,'r','linewidth',1.5);
title('Fixed attenuation channel')

xlabel('Transmissivity $T_E$')
ylabel('Rate (bit/pulse)')
legend('No-PS', 'T-PS', 'R-PS', 'location', 'southeast')
ylim([1e-10, 5])
xlim([0, 1])

axes('position',[0.25,0.2,0.4,0.3]);
semilogy(10.^((DistanceNo - Attenuation)/-50), RateNo ,'k--','linewidth',1.5);
hold on
semilogy(10.^((DistanceT - Attenuation)/-50),RateT,'b','linewidth',1.5);
semilogy(10.^((DistanceR - Attenuation)/-50),RateR,'r','linewidth',1.5);
axis([0.005,0.02,1e-10,0.005])
saveas(gcf, 'Fig3', 'epsc')

figure(2)
semilogy(DistanceNo - Attenuation, RateNo ,'k--','linewidth',1.5)
hold on
semilogy(DistanceT' - Attenuation, RateT','b','linewidth',1.5);
semilogy(DistanceR' - Attenuation, RateR','r','linewidth',1.5);
title('Fixed attenuation channel')

xlabel('Distance($km$)')
ylabel('Rate (bit/pulse)')
legend('No-PS', 'T-PS', 'R-PS', 'location', 'southwest')
ylim([1e-10, 5])
xlim([0, 120])
saveas(gcf, 'Fig4', 'epsc')

% %% 3D plot
% % Plot Rate vs. alpha^2(A)
% ParaPho.mode = 'CalRat';
% figure(3)
% Distance = linspace(60,115,30);
% T_Range = 10.^(-0.2*Distance/10);
% 
% zP = 1;
% B = 0.001;
% Abin = linspace(0.05,5,35);
% Alen = length(Abin);
% 
% NoSubtraction_Data = zeros(Alen, length(Distance));
% RxSubtraction_Data = NoSubtraction_Data;
% TxSubtraction_Data = NoSubtraction_Data;
% 
% tic
% parfor a = 1:Alen    
%     NoSubtraction_Data(a,:) = Operation(T_Range, 1,...
%         struct('B', B, 'th', th, 'Tbs', Tbs, 'sum_limit_n', sum_limit_n,...
%         'sum_limit_m', sum_limit_m, 'mode', 'CalRat', 'operation', 'No'),...
%         Abin(a), z, 0);
%     
%     RxSubtraction_Data(a,:) = Operation(T_Range, 1,...
%         struct('B', B, 'th', th, 'Tbs', Tbs, 'sum_limit_n', sum_limit_n,...
%         'sum_limit_m', sum_limit_m, 'mode', 'CalRat', 'operation', 'RxS'),...
%         Abin(a), z, 0);
%     
%     TxSubtraction_Data(a,:) = Operation(T_Range, 1,...
%         struct('B', B, 'th', th, 'Tbs', Tbs, 'sum_limit_n', sum_limit_n,...
%         'sum_limit_m', sum_limit_m, 'mode', 'CalRat', 'operation', 'TxS'),...
%         Abin(a), z, 0);
% end
% toc
% 
% [~,minIndexR] = min(RxSubtraction_Data,[],2);
% [~,minIndexT] = min(TxSubtraction_Data,[],2);
% [~,minIndexN] = min(NoSubtraction_Data,[],2);
% for a = 1:Alen
%     if minIndexR(a) < length(Distance) - zP
%         RxSubtraction_Data(a,minIndexR(a)+zP:end) = nan;
%     end
%     if minIndexT(a) < length(Distance) - zP
%         TxSubtraction_Data(a,minIndexT(a)+zP:end) = nan;
%     end
%     if minIndexN(a) < length(Distance) - zP
%         NoSubtraction_Data(a,minIndexN(a)+zP:end) = nan;
%     end
% end
% 
% surf(Distance, Abin, RxSubtraction_Data,'Edgecolor','r','linewidth',.5)
% hold on
% surf(Distance, Abin, TxSubtraction_Data,'Edgecolor','b','linewidth',.5)
% surf(Distance, Abin, NoSubtraction_Data,'Edgecolor','k','linewidth',.5)
% surf([min(Distance) max(Distance)],[min(Abin) max(Abin)],zeros(2,2),'Edgecolor','n')
% xlim([min(Distance) max(Distance)])
% ylim([min(Abin) max(Abin)])
% grid on
% box on
% colormap('cool')
% xlabel('Distance ($km$)')
% ylabel('$\alpha^2$')
% zlabel('Rate (bit/pulse)')
% title('Rate vs. mean photon number $\alpha^2$')
% view([-137,31])
% % legend('R-PS', 'T-PS', 'No-PS','location','best')
% saveas(gcf, 'Fig6', 'epsc')

% %%
% figure(4)
% Distance = linspace(50,100,30);
% T_Range = 10.^(-0.2*Distance/10);
% 
% A = 1.3;
% Bbin = linspace(0.001,0.01,30);
% Blen = length(Bbin);
% 
% NoSubtraction_Data = zeros(Blen, length(Distance));
% RxSubtraction_Data = NoSubtraction_Data;
% TxSubtraction_Data = NoSubtraction_Data;
% Ground = NoSubtraction_Data;
% tic
% parfor b = 1:Blen    
%     NoSubtraction_Data(b,:) = Operation(T_Range, 1,...
%         struct('B', Bbin(b), 'th', th, 'Tbs', Tbs, 'sum_limit_n', sum_limit_n,...
%         'sum_limit_m', sum_limit_m, 'mode', 'CalRat', 'operation', 'No'),...
%         A, z, 0);
%     
%     RxSubtraction_Data(b,:) = Operation(T_Range, 1,...
%         struct('B', Bbin(b), 'th', th, 'Tbs', Tbs, 'sum_limit_n', sum_limit_n,...
%         'sum_limit_m', sum_limit_m, 'mode', 'CalRat', 'operation', 'RxS'),...
%         A, z, 0);
%     
%     TxSubtraction_Data(b,:) = Operation(T_Range, 1,...
%         struct('B', Bbin(b), 'th', th, 'Tbs', Tbs, 'sum_limit_n', sum_limit_n,...
%         'sum_limit_m', sum_limit_m, 'mode', 'CalRat', 'operation', 'TxS'),...
%         A, z, 0);
% end
% toc
% 
% zP = 1;
% [~,minIndexR] = min(RxSubtraction_Data,[],2);
% [~,minIndexT] = min(TxSubtraction_Data,[],2);
% [~,minIndexN] = min(NoSubtraction_Data,[],2);
% for a = 1:Blen
%     if minIndexR(a) < length(Distance) - zP
%         RxSubtraction_Data(a,minIndexR(a)+zP:end) = nan;
%     end
%     if minIndexT(a) < length(Distance) - zP
%         TxSubtraction_Data(a,minIndexT(a)+zP:end) = nan;
%     end
%     if minIndexN(a) < length(Distance) - zP
%         NoSubtraction_Data(a,minIndexN(a)+zP:end) = nan;
%     end
% end
% 
% surf(Distance, Bbin, RxSubtraction_Data,'Edgecolor','r','linewidth',.5)
% hold on
% surf(Distance, Bbin, TxSubtraction_Data,'Edgecolor','b','linewidth',.5)
% surf(Distance, Bbin, NoSubtraction_Data,'Edgecolor','k','linewidth',.5)
% surf([min(Distance) max(Distance)],[min(Bbin) max(Bbin)],zeros(2,2),'Edgecolor','n')
% xlim([min(Distance) max(Distance)])
% ylim([min(Bbin) max(Bbin)])
% grid on
% box on
% colormap('cool')
% xlabel('Distance (km)')
% ylabel('$\beta^2$')
% zlabel('Rate (bit/pulse)')
% title('Rate vs. channel noise $\beta^2$')
% % legend('R-PS', 'T-PS', 'No-PS','location','best')
% view(-133,37)
% saveas(gcf, 'Fig5', 'epsc')

% %% Avg Rate
% figure(5)
% ParaPho.mode = 'CalRat';
% sigmaBin = [linspace(0.001,1,20) linspace(1.0001,20,20)];
% sigmaLen = length(sigmaBin);
% 
% Av_Rate_No = zeros(length(sigmaBin), 1);
% Av_Rate_TxS = Av_Rate_No;
% Av_Rate_RxS = Av_Rate_No;
% 
% eta0 = sqrt(1-exp(-2)); % upper limit of integral
% 
% tic
% ParaPho.operation = 'No';
% parfor sg = 1:sigmaLen
%     %disp(s/length(sigmaBin)) % progress indicator
%     Av_Rate_No(sg) = integral(@(x) WeibullPDF(x,sigmaBin(sg)) .* Operation(x.^2, 1, ParaPho, A, z, 0),0.001,eta0);    
% end
% 
% ParaPho.operation = 'TxS';
% parfor sg = 1:sigmaLen
%     Av_Rate_TxS(sg) = integral(@(x) WeibullPDF(x,sigmaBin(sg)) .* Operation(x.^2, 1, ParaPho, A, z, 0),0.001,eta0);    
% end
% 
% ParaPho.operation = 'RxS';
% parfor sg = 1:sigmaLen
%     Av_Rate_RxS(sg) = integral(@(x) WeibullPDF(x,sigmaBin(sg)) .* Operation(x.^2, 1, ParaPho, A, z, 0),0.001,eta0);
% end
% toc
% 
% semilogy(sigmaBin, Av_Rate_No,'k--',...
%     sigmaBin, Av_Rate_TxS,'b',sigmaBin, Av_Rate_RxS,'r',sigmaBin,...
%     Av_Rate_TxS*Av_Rate_No(1)/Av_Rate_TxS(1),'b--','LineWidth', 1)
% legend('No-PS','T-PS','R-PS', 'T-PS Normalized','Location', 'northeast')
% title(['Average rate when \alpha^2=' num2str(A) ' and \beta^2=' num2str(B)])
% ylabel('Average rate (bit/pulse)')
% xlabel('\sigma_b')
% saveas(gcf, 'Fig7', 'epsc')

% %%
% figure(6)
% sigmaBin = sigmaBin(1:20);
% sigmaLen = length(sigmaBin);
% semilogy(sigmaBin, Av_Rate_No(1:20),'k--',...
%     sigmaBin, Av_Rate_TxS(1:20),'b',sigmaBin, Av_Rate_RxS(1:20),'r',sigmaBin,...
%     Av_Rate_TxS(1:20)*Av_Rate_No(1)/Av_Rate_TxS(1),'b--','LineWidth', 1)
% legend('No-PS','T-PS','R-PS', 'T-PS Normalized','Location', 'southwest')
% hold on
% ylabel('Average rate (bit/pulse)')
% xlabel('\sigma_b')
% 
% Av_Rate_No = zeros(length(sigmaBin), 1);
% Av_Rate_TxS = Av_Rate_No;
% Av_Rate_RxS = Av_Rate_No;
% 
% tic
% ParaPho.B = .5;
% ParaPho.operation = 'No';
% parfor sg = 1:sigmaLen
%     %disp(s/length(sigmaBin)) % progress indicator
%     Av_Rate_No(sg) = integral(@(x) WeibullPDF(x,sigmaBin(sg)) .* Operation(x.^2, 1, ParaPho, A, z, 0),0.001,eta0);    
% end
% 
% ParaPho.operation = 'TxS';
% parfor sg = 1:sigmaLen
%     Av_Rate_TxS(sg) = integral(@(x) WeibullPDF(x,sigmaBin(sg)) .* Operation(x.^2, 1, ParaPho, A, z, 0),0.001,eta0);    
% end
% 
% ParaPho.operation = 'RxS';
% parfor sg = 1:sigmaLen
%     Av_Rate_RxS(sg) = integral(@(x) WeibullPDF(x,sigmaBin(sg)) .* Operation(x.^2, 1, ParaPho, A, z, 0),0.001,eta0);
% end
% toc
% semilogy(sigmaBin, Av_Rate_No(1:20),'k--',...
%     sigmaBin, Av_Rate_TxS(1:20),'b',sigmaBin, Av_Rate_RxS(1:20),'r',sigmaBin,...
%     Av_Rate_TxS(1:20)*Av_Rate_No(1)/Av_Rate_TxS(1),'b--','LineWidth', .5)
% legend('No-PS','T-PS','R-PS', 'T-PS Normalized','Location', 'southwest')
% title(['Average rate when \alpha^2=' num2str(A)])
% saveas(gcf, 'Fig8', 'epsc')
