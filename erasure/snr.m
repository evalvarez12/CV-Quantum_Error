function [res] = snr(V, sigma, T1r, T2r, T3r, D1, D2)
    T1 = sqrt(T1r);
    T2 = sqrt(T2r);
    T3 = sqrt(T3r);

    a = V/4*((T1+T2)^2/8 + T3^2/2) - sqrt(V^1-1)*(T1+T2)*T3/8; 
    b = (T1-T2)^2/32;
    c = (1 - T1/4 - T2/4 - T3/2)/4;
    d = (T1-T2)/2;
    ap = a+b+c;

    sigma_snr = sqrt(d^2*sigma + 4*ap);
    mu_snr = (T2+T1)*D1/2 + T3*D2;
    res = sigma_snr/mu_snr;
end