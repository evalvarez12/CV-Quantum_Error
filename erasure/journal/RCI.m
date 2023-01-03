function [res]=RCI(V,vs, Tr, eps, eta, pe)
    % Same T for all the modes
%     gp = g*eta;

%     T1 = sqrt(T1r);
%     T2 = sqrt(T2r);
%     T3 = sqrt(T3r);
% 
%     Tm = (T1 - T2)/2;
%     Tp = (T1 + T2)/2;

%     e1 = V*((Tm+gp*Tp)^2 + (T3*gp)^2) - 2*V2*(Tm+gp*Tp)*T3*gp;
%     e2 = 1 + g.^2*(1+eta^2) + eps*(1+2.*gp.^2) - T1^2/2.*(1+gp).^2 - T2^2/2.*(1-gp).^2 - T3^2.*gp.^2;
% 
%     e10 = 0;
%     e20 = 1 + eps - Tr 

    T = sqrt(Tr);
    V2 = sqrt(V^2 - 1);
    % No erasure CM
    a1 = vs;
    b1 = Tr*(vs-1) + eps + 1;
    c1 = T*sqrt(vs^2-1);

    z = sqrt((a1+b1).^2 - 4*c1.^2);
    vne1 = (z+(a1-b1))/2;
    vne2 = (z-(a1-b1))/2;

    % Erasure CM --- g=1/eta
    a2 = vs;    
    b2 = Tr*vs + 2*Tr*(V-V2) + 3 + 3*eps - 3*Tr + (1/eta^2-1);
    c2 = T*sqrt(vs^2-1);

    z2 = sqrt((a2+b2).^2 - 4*c2.^2);
    ve1 = (z2+(a2-b2))/2;
    ve2 = (z2-(a2-b2))/2;
    

    Sth = gfun(vs);
    Sne = real(gfun(vne1) + gfun(vne2));
    Se = real(gfun(ve1) + gfun(ve2));

    Se(Se<0) = 0;
    Sne(Sne<0) = 0;
    
    Sth
    Se
    Sne

    VV = V-V2;
    VV

    res= (1 - 3*pe^2*(1-pe) - pe^3)*Sth - ((1-pe)^3 + pe*(1-pe)^2)*Sne - 2*pe*(1-pe)^2*Se;
end