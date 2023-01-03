clear all;

vs = 5e12;
T=0.8;
eps = 0.05;
V= 10;
V2 = sqrt(V^2-1);
eta =1;


a = vs
b = T.*(vs-1) + eps +1
c = sqrt(T).*sqrt(vs^2-1)

z = sqrt((a+b).^2 - 4*c.^2)
v1 = (z+(a-b))/2
v2 = (z-(a-b))/2


a2 = vs;    
b2 = T*vs + 2*T*(V-V2) + 3 + 3*eps - 3*T + (1/eta^2-1);
c2 = sqrt(T)*sqrt(vs^2-1);

z2 = sqrt((a2+b2).^2 - 4*c2.^2);
ve1 = (z2+(a2-b2))/2;
ve2 = (z2-(a2-b2))/2;
    



Sth = gfun(vs)
S = real(gfun(v1) + gfun(v2))
Se = real(gfun(ve1) + gfun(ve2))

Sth=1;
S=1;
Se=1;

pe=0.6;
PRO =(1 - 3*pe^2*(1-pe) - pe^3)*Sth - ((1-pe)^3 + pe*(1-pe)^2)*S - 2*pe*(1-pe)^2*Se
DIR= (1-pe)*Sth  - (1-pe)*S
