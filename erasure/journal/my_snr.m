function [snr] = my_snr(T1, T2, T3, sigma, V, delta)
    V2 = sqrt(V^2 -1);
    
    a = V*((T1 + T2)^2/8 + T3^2/2) - 2*V2*(T1 + T2)*T3/4;
    b = (T1 - T2)^2/8;
    c = (1 - T1^2/4 - T2^2/4 - T3^2/2);
    d = (T1 - T2)/2;
a
b
c
d
    epsilon = 0.1;
    sig_snr = d^2*sigma + (a + b + c) + epsilon;
%     snr = (T3 * delta)/ (sig_snr); 
    mu_snr =  ((T1+T2)/2 * delta);
sig_snr
mu_snr
    snr = (mu_snr)/ (sig_snr); 
%     snr = (((T1+T2)/2 + T3)/2 * delta)/ (sig_snr); 
end