clear all;
close all;

% PD
X = 0:0.005:0.5;

data = load('F_pd_mean.mat');
Fm = data.Fm;

data = load('F_pd_std.mat');
Fstd = data.Fstd;

data = load('F_pd_mean_dir.mat');
Fm_dir = data.Fm_dir;

data = load('F_pd_std_dir.mat');
Fstd_dir = data.Fstd_dir;

figure;
hold on;
errorbar(X, Fm, Fstd);
errorbar(X, Fm_dir, Fstd_dir);

% T0
X = 0:0.02:0.8;

data = load('F_T0_mean.mat');
Fm = data.Fm;

data = load('F_T0_std.mat');
Fstd = data.Fstd;

data = load('F_T0_mean_dir.mat');
Fm_dir = data.Fm_dir;

data = load('F_T0_std_dir.mat');
Fstd_dir = data.Fstd_dir;

figure;
hold on;
errorbar(X, Fm, Fstd);
errorbar(X, Fm_dir, Fstd_dir);

% sigma
X = 0:0.002:0.2;

data = load('F_sigma_mean.mat');
Fm = data.Fm;

data = load('F_sigma_std.mat');
Fstd = data.Fstd;

data = load('F_sigma_mean_dir.mat');
Fm_dir = data.Fm_dir;

data = load('F_sigma_std_dir.mat');
Fstd_dir = data.Fstd_dir;

figure;
hold on;
errorbar(X, Fm, Fstd);
errorbar(X, Fm_dir, Fstd_dir);