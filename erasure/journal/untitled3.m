clear all;

vs = 400;
T=1;
eps = 0.0;
V= 50;
V2 = sqrt(V^2-1);



ai = vs;
bi = vs;
ci = sqrt(vs^2-1);


adir = vs;
bdir = T.*(vs-1) + eps +1;
cdir = sqrt(T).*sqrt(vs^2-1);

a = vs;    
b = T*vs + 2*T*(V-V2) + 3 + 3*eps - 3*T;
c = sqrt(T)*sqrt(vs^2-1);

f_tmsv_gauss(ai, bi, ci, a, b, c)

% f_code_tmsv(0,1,1,0,-1,10,10)
% 
% ans =
% 
%     0.6662




% Sth = gfun(vs)
% S = real(gfun(v1) + gfun(v2))
% Se = real(gfun(ve1) + gfun(ve2))
% 
% Sth=1;
% S=1;
% Se=1;
% 
% pe=0.6;
% PRO =(1 - 3*pe^2*(1-pe) - pe^3)*Sth - ((1-pe)^3 + pe*(1-pe)^2)*S - 2*pe*(1-pe)^2*Se
% DIR= (1-pe)*Sth  - (1-pe)*S
