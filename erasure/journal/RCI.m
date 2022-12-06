function [res]=RCI(V, g,vs, T1r, T2r, T3r, eps, eta)

    gp = g*eta;

    T1 = sqrt(T1r);
    T2 = sqrt(T2r);
    T3 = sqrt(T3r);

    Tm = (T1 - T2)/2;
    Tp = (T1 + T2)/2;
    V2 = sqrt(V^2 - 1);

    e1 = V*((Tm+gp*Tp)^2 + (T3*gp)^2) - 2*V2*(Tm+gp*Tp)*T3*gp;
    e2 = 1 + g.^2*(1+eta^2) + eps*(1+2.*gp.^2) - T1^2/2.*(1+gp).^2 - T2^2/2.*(1-gp).^2 - T3^2.*gp.^2;

    
    a = vs;
    b = (Tp + gp*Tm)^2*vs + e1 + e2;
    c = (Tp + gp*Tm)*sqrt(vs^2-1);

    z = sqrt((a+b)^2 - 4*c^2);
    v1 = (z+(b-a))/2;
    v2 = (z-(b-a))/2;
    v3 = vs;

    res=gfun(v3) - gfun(v1) - gfun(v2);

end