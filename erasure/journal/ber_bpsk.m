function [ber] = ber_bpsk(T1, T2, T3, sigma, V, delta)
    [mu_snr, sig_snr] = my_snr(T1, T2, T3, sigma, V, delta);

    Q = mu_snr / (sig_snr);
    
    ber = 1/2 * erfc(Q/sqrt(2));
 
end