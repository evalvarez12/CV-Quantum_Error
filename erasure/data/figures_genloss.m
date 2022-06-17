clear all;
close all;

% PD
x = 0:0.005:0.5;
X = repmat(x,200,1);
X = transpose(X);

data = load('F_pd.mat');
Fm = data.Fm;

% data = load('F_pd_std.mat');
% Fstd = data.Fstd;

data = load('F_pd_dir.mat');
Fm_dir = data.Fm_dir;

% data = load('F_pd_std_dir.mat');
% Fstd_dir = data.Fstd_dir;

figure;
hold on;
scatter(X, Fm, 1, 'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1);

figure;
hold on;
scatter(X, Fm_dir, 1, '+' ,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1);
% errorbar(X, Fm, Fstd);
% errorbar(X, Fm_dir, Fstd_dir);

% % T0
% X = 0:0.02:0.8;
% 
% data = load('F_T0_mean.mat');
% Fm = data.Fm;
% 
% data = load('F_T0_std.mat');
% Fstd = data.Fstd;
% 
% data = load('F_T0_mean_dir.mat');
% Fm_dir = data.Fm_dir;
% 
% data = load('F_T0_std_dir.mat');
% Fstd_dir = data.Fstd_dir;
% 
% figure;
% hold on;
% errorbar(X, Fm, Fstd);
% errorbar(X, Fm_dir, Fstd_dir);
% 
% % sigma
% X = 0:0.002:0.2;
% 
% data = load('F_sigma_mean.mat');
% Fm = data.Fm;
% 
% data = load('F_sigma_std.mat');
% Fstd = data.Fstd;
% 
% data = load('F_sigma_mean_dir.mat');
% Fm_dir = data.Fm_dir;
% 
% data = load('F_sigma_std_dir.mat');
% Fstd_dir = data.Fstd_dir;
% 
% figure;
% hold on;
% errorbar(X, Fm, Fstd);
% errorbar(X, Fm_dir, Fstd_dir);