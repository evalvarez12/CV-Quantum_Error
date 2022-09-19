function [snr] = SNR(T1, T2, T3, sigma, V, delta)
    V2 = sqrt(V^2 -1);
    
    a = 1/4*V*((T1 + T2)^2/8 + T3^2/2) - V2/2*(T1 + T2)*T3/4;
    b = 1/4*(T1 - T2)^2/8;
    c = 1/4*(1 - T1^2/4 - T2^2/4 - T3^2/2);
    d = (T1 - T2)/2;

    epsilon = 0.1;
    sig_snr = d^2*sigma + 4*(a + b + c) + epsilon;
%     snr = (T3 * delta)/ (sig_snr); 
    snr = ((T1+T2)/2 * delta)/ (sig_snr); 
%     snr = (((T1+T2)/2 + T3)/2 * delta)/ (sig_snr); 
end