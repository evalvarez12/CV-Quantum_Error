% file name : WeibullPDF.m
function [ T,eta ] = WeibullPDF( eta, sigma_b, beta, W )
%% Do not use sigma_b = 0
%% Default parameters
if nargin < 3
    beta = 1;
    W = 1;
    h = (beta / W)^2;
%     eta0 = sqrt (1 - exp(-2 * h));
%     eta = 0:0.00001:(eta0-0.00001);
%     sigma_b = 0.5;
end
h = (beta / W)^2;
eta0 = sqrt(1 - exp(-2 * h));
% PDF related parameter

% I0 = besseli(0, 4 * h);
% I1 = besseli(1, 4 * h);
% Exph = exp(-4 * h);
% Rh = log(2 * eta0^2 / (1 - Exph * I0));
% lambda = 8 * h * Exph * I1 / (1 - Exph * I0) * Rh^(-1);
% L = beta * Rh^(-1/lambda);

lambda = 2.3129;
L = 1.1136;

%% Weibull PDF
if sigma_b
T = 2 * L^2 / sigma_b^2 / lambda ./eta ...
         .* (2 * log(eta0./eta)).^(2/lambda - 1)...
         .* exp(-L^2 / 2 / sigma_b^2 * ((2 * log(eta0./eta)).^(2/lambda)) );
    T(eta==0) = 0;
else
    T = zeros(1,length(eta));
end

