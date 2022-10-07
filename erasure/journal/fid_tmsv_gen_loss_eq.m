function [F] = fid_tmsv_gen_loss_eq(V, g, sigma, T1r, T2r, T3r, eps, eta)


gp = g*eta;

T1 = sqrt(T1r);
T2 = sqrt(T2r);
T3 = sqrt(T3r);

T12 = (T1 - T2)/2;
T21 = (T1 + T2)/2;
V2 = sqrt(V^2 - 1);

T1p = sqrt(1 + eps - T1r);
T2p = sqrt(1 + eps - T2r);
T3p = sqrt(1 + eps - T3r);


a = (V*((T12 + gp.*T21).^2 + (T3.*gp).^2)) - 2*V2*(T12 + gp.*T21).*T3.*gp;
b = (T21 + gp.*T12).^2;
c = 1 + g.^2*(1+eta^2) + eps*(1+2.*gp.^2) - T1^2/2.*(1+gp).^2 - T2^2/2.*(1-gp).^2 - T3^2.*gp.^2;

d = sqrt(2)*(T21 + gp.*T12);

a2 = a + b + c + 1;
b2 = d - sqrt(2);

% F =  4 ./ (b2.^2*sigma + 4*a2);
F =  2 ./ (b2.^2*sigma + a2);

% F = real(F);
end