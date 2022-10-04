function [F] = fid_tmsv_gen_loss_eq_old(V, g, sigma, T1r, T2r, T3r, eps)

T1 = sqrt(T1r);
T2 = sqrt(T2r);
T3 = sqrt(T3r);

T12 = (T1 - T2)/2;
T21 = (T1 + T2)/2;
V2 = sqrt(V^2 - 1);

T1p = sqrt(1 + eps - T1r);
T2p = sqrt(1 + eps - T2r);
T3p = sqrt(1 + eps - T3r);


a = 1/4*(V*((T12 + g.*T21).^2 + (T3.*g).^2)) - V2/2*(T12 + g.*T21).*T3.*g;
b = 1/4*(T21 + g.*T12).^2;
c = 1/4*(T1p.^2/2 + T2p.^2/2 + g.*(T1p.^2 - T2p.^2) + g.^2.*(T1p.^2/2 + T2p.^2/2 + T3p.^2));
d = sqrt(2)*(T21 + g.*T12);

a2 = a + b + c + 1/4;
b2 = d - sqrt(2);

% F =  4 ./ (b2.^2*sigma + 4*a2);
F =  2 ./ (b2.^2*sigma + 4*a2);

% F = real(F);
end