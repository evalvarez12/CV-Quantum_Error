clear all;
close all;
set(0,'defaultTextInterpreter','latex');

% Make means and save
% data1 = load('Fscatter_V3_sigma2.mat');
% data1 = data1.Fm;
% Fm = mean(data1,2);
% save('Fscatter_mean_V3_sigma2.mat', "Fm")
% 
% 
data2 = load('Ferasure_T1_V6_sigma2.mat');
data2 = data2.Fm;
Fm = mean(data2,2);
save('Ferasure_mean_T1_V6_sigma2.mat', "Fm")


par_erasure = linspace(.5,1,24);
par_no = .5:0.002:1;



data1_erasure = load('Ferasure_mean_T1_V3_sigma2.mat');
data1_erasure = transpose(data1_erasure.Fm);

data2_erasure = load('Ferasure_mean_T1_V6_sigma2.mat');
data2_erasure = transpose(data2_erasure.Fm);


data3_erasure = load('Ferasure_mean_T1_V10_sigma2.mat');
data3_erasure = transpose(data3_erasure.Fm);

data1_no = load('Fscatter_mean_V3_sigma2.mat');
data1_no = transpose(data1_no.Fm);

data2_no = load('Fscatter_mean_V10_sigma2.mat');
data2_no = transpose(data2_no.Fm);



sigma_coh = 10;
F_class = (1 + (1/sigma_coh))/(2 + (1/sigma_coh));



set(gca,'fontname','times') 

figure;
hold on;
plot(par_erasure, data3_erasure, 'p-', 'LineWidth', 1.5, 'DisplayName', 'V=10');
plot(par_erasure, data2_erasure, '*-', 'LineWidth', 1.5, 'DisplayName', 'V=6');
plot(par_erasure, data1_erasure, 'v-', 'LineWidth', 1.5, 'DisplayName', 'V=3');


% plot(par_no, data1_no, 'LineWidth', 1.2, 'DisplayName', 'V=3');
plot(par_no, data2_no, 'LineWidth', 1.5, 'DisplayName', 'No erasure');

ylim([0.4 1]);

xlabel('$T$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
% legend('Location','northwest')
text(0.95, 0.97, '$\sigma=0.2$')

savefigures('erasure_sigma2')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data1_erasure = load('Ferasure_T1_V3_sigma0.mat');
data1_erasure = transpose(data1_erasure.Fm);

data2_erasure = load('Ferasure_T1_V6_sigma0.mat');
data2_erasure = transpose(data2_erasure.Fm);


data3_erasure = load('Ferasure_T1_V10_sigma0.mat');
data3_erasure = transpose(data3_erasure.Fm);

data1_no = load('Fnoerasure_V3_sigma0.mat');
data1_no = transpose(data1_no.Fm);




figure;
hold on;
plot(par_erasure, data3_erasure, 'p-', 'LineWidth', 1.5, 'DisplayName', 'V=10');
plot(par_erasure, data2_erasure, '*-', 'LineWidth', 1.5, 'DisplayName', 'V=6');
plot(par_erasure, data1_erasure, 'v-', 'LineWidth', 1.5, 'DisplayName', 'V=3');


% plot(par_no, data1_no, 'LineWidth', 1.2, 'DisplayName', 'V=3');
plot(par_erasure, data1_no, 'LineWidth', 1.5, 'DisplayName', 'No erasure');



xlabel('$T$', 'Interpreter', 'latex');
ylabel('$\bar{\mathcal{F}}$', 'Interpreter', 'latex');
legend('Location','northwest')
text(0.81, 0.97, '$\sigma=0$')


savefigures('erasure_sigma0')





