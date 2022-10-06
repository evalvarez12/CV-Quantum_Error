function [mu2_snr, sig2_snr] = my_snr(T1r, T2r, T3r, sigma, V, eps, etar, delta)
    V2 = sqrt(V^2 -1);

    T1 = sqrt(T1r);
    T2 = sqrt(T2r);
    T3 = sqrt(T3r);
    eta = sqrt(etar);


    Tp = (T1 + T2)/2;
    Tm = (T1 - T2)/2;
    
    a = etar*V*(Tp^2/2 + T3^2/2) - 2*V2*Tp*T3/2;
    b = etar*Tm^2/2;
    c = 1 + etar*(eps - T1^2/4 - T2^2/4 - T3^2/2);
    d = eta*Tm;

    sig2_snr = d^2*sigma + (a + b + c);
%     snr = (T3 * delta)/ (sig_snr); 
    mu2_snr =  (etar*Tp^2 * delta^2 );

%     snr = (mu_snr)/ (sig_snr); 
%     snr = (((T1+T2)/2 + T3)/2 * delta)/ (sig_snr); 
end