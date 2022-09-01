function [F] = f_code(T1, T2, T3, V, T, g, operation, sigma)
warning('off')

switch operation
    case 'tmsv'
        F=...
            integralN(@(r,theta,ra,ri)...
            exp(-(ra.^2 + ri.^2)./sigma).*1/(pi*sigma).*...
            (exp(1).^((-1/8).*r.*(4.*r+8.*exp(1).^(sqrt(-1).*theta).*(ra+(sqrt( ...
  -1)*(-1)).*ri)+(-8).*exp(1).^((sqrt(-1)*(-1)).*theta).*(ra+sqrt( ...
  -1).*ri)+(-4).*exp(1).^(sqrt(-1).*theta).*(ra+(sqrt(-1)*(-1)).*ri) ...
  .*(T1+g.*T1+T2+(-1).*g.*T2)+4.*exp(1).^((sqrt(-1)*(-1)).*theta).*( ...
  ra+sqrt(-1).*ri).*(T1+g.*T1+T2+(-1).*g.*T2)+r.*(T1+g.*T1+T2+(-1).* ...
  g.*T2).^2+(-2).*r.*((-2)+T1.^2+T2.^2+2.*g.*(T1.^2+(-1).*T2.^2)+ ...
  g.^2.*((-4)+T1.^2+T2.^2+2.*T3.^2))+4.*r.*((1/4).*((1+g).*T1+((-1)+ ...
  g).*T2).^2+g.^2.*T3.^2).*V+(-4).*g.*r.*((1+g).*T1+((-1)+g).*T2).* ...
  T3.*((-1)+V.^2).^(1/2))).*pi.^(-1).*r)...
            ,0,inf,0,2*pi,-3*sigma,3*sigma,-3*sigma,3*sigma);
    
    case 'ps'
        F=...
            integralN(@(r,theta,ra,ri)...
            exp(-(ra.^2 + ri.^2)./sigma).*1/(pi*sigma).*...
            ((1/4).*exp(1).^((-1/2).*r.*(r+2.*exp(1).^(sqrt(-1).*theta).*(ra+( ...
  sqrt(-1)*(-1)).*ri)+(-2).*exp(1).^((sqrt(-1)*(-1)).*theta).*(ra+ ...
  sqrt(-1).*ri)+(-1).*exp(1).^(sqrt(-1).*theta).*(ra+(sqrt(-1)*(-1)) ...
  .*ri).*(T1+g.*T1+T2+(-1).*g.*T2)+exp(1).^((sqrt(-1)*(-1)).*theta) ...
  .*(ra+sqrt(-1).*ri).*(T1+g.*T1+T2+(-1).*g.*T2)+(1/4).*r.*(T1+g.* ...
  T1+T2+(-1).*g.*T2).^2+(-1/2).*r.*((-2)+T1.^2+T2.^2+2.*g.*(T1.^2+( ...
  -1).*T2.^2)+g.^2.*((-4)+T1.^2+T2.^2+2.*T3.^2))+(-1).*g.^2.*r.* ...
 T3.^2.*((-1)+T.*((-1)+V)+(-1).*V).^(-1).*(1+T.*((-1)+V)+V)+2.*g.* ...
  r.*((1+g).*T1+((-1)+g).*T2).*T3.*((-1)+T.*((-1)+V)+(-1).*V).^(-1) ...
  .*(T.*((-1)+V.^2)).^(1/2)+(-1).*r.*((1+g).*T1+((-1)+g).*T2).^2.*( ...
  1+T.*((-1)+V)+V).*(4.*T.*((-1)+V)+(-4).*(1+V)).^(-1))).*pi.^(-1).* ...
  r.*((-1)+T.*((-1)+V)+(-1).*V).^(-1).*(4.*T.*(1+g.^2.*r.^2.*T3.^2) ...
  .*((-1)+V)+(-4).*(1+V)+r.^2.*((1+g).*T1+((-1)+g).*T2).*((1+g).* ...
  T1.*(1+V)+((-1)+g).*T2.*(1+V)+(-4).*g.*T3.*(T.*((-1)+V.^2)).^(1/2) ...
  )))...
            ,0,inf,0,2*pi,-3*sigma,3*sigma,-3*sigma,3*sigma);
    
    case 'pa'
        F=...
            integralN(@(r,theta,ra,ri)...
            exp(-(ra.^2 + ri.^2)./sigma).*1/(pi*sigma).*...
            ((1/4).*exp(1).^((-1/2).*r.*(r+2.*exp(1).^(sqrt(-1).*theta).*(ra+( ...
  sqrt(-1)*(-1)).*ri)+(-2).*exp(1).^((sqrt(-1)*(-1)).*theta).*(ra+ ...
  sqrt(-1).*ri)+(-1).*exp(1).^(sqrt(-1).*theta).*(ra+(sqrt(-1)*(-1)) ...
  .*ri).*(T1+g.*T1+T2+(-1).*g.*T2)+exp(1).^((sqrt(-1)*(-1)).*theta) ...
  .*(ra+sqrt(-1).*ri).*(T1+g.*T1+T2+(-1).*g.*T2)+(1/4).*r.*(T1+g.* ...
  T1+T2+(-1).*g.*T2).^2+(-1/2).*r.*((-2)+T1.^2+T2.^2+2.*g.*(T1.^2+( ...
  -1).*T2.^2)+g.^2.*((-4)+T1.^2+T2.^2+2.*T3.^2))+(-1).*g.^2.*r.* ...
 T3.^2.*((-1)+T.*((-1)+V)+(-1).*V).^(-1).*(1+T.*((-1)+V)+V)+2.*g.* ...
  r.*((1+g).*T1+((-1)+g).*T2).*T3.*((-1)+T.*((-1)+V)+(-1).*V).^(-1) ...
  .*(T.*((-1)+V.^2)).^(1/2)+(-1).*r.*((1+g).*T1+((-1)+g).*T2).^2.*( ...
  1+T.*((-1)+V)+V).*(4.*T.*((-1)+V)+(-4).*(1+V)).^(-1))).*pi.^(-1).* ...
  r.*((-1)+T.*((-1)+V)+(-1).*V).^(-1).*((-4)+T.*(4+r.^2.*((1+g).*T1+ ...
  ((-1)+g).*T2).^2).*((-1)+V)+(-4).*V+(-4).*g.*r.^2.*(T1+(-1).*T2).* ...
  T3.*(T.*((-1)+V.^2)).^(1/2)+4.*g.^2.*r.^2.*T3.*(T3.*(1+V)+(-1).*( ...
  T1+T2).*(T.*((-1)+V.^2)).^(1/2))))...
            ,0,inf,0,2*pi,-3*sigma,3*sigma,-3*sigma,3*sigma);
    
    case 'pc'
        F=...
            integralN(@(r,theta,ra,ri)...
            exp(-(ra.^2 + ri.^2)./sigma).*1/(pi*sigma).*...
            (exp(1).^((-1/2).*r.*(r+2.*exp(1).^(sqrt(-1).*theta).*(ra+(sqrt(-1) ...
  *(-1)).*ri)+(-2).*exp(1).^((sqrt(-1)*(-1)).*theta).*(ra+sqrt(-1).* ...
  ri)+(-1).*exp(1).^(sqrt(-1).*theta).*(ra+(sqrt(-1)*(-1)).*ri).*( ...
  T1+g.*T1+T2+(-1).*g.*T2)+exp(1).^((sqrt(-1)*(-1)).*theta).*(ra+ ...
  sqrt(-1).*ri).*(T1+g.*T1+T2+(-1).*g.*T2)+(1/4).*r.*(T1+g.*T1+T2+( ...
  -1).*g.*T2).^2+(-1/2).*r.*((-2)+T1.^2+T2.^2+2.*g.*(T1.^2+(-1).* ...
  T2.^2)+g.^2.*((-4)+T1.^2+T2.^2+2.*T3.^2))+(-1).*g.^2.*r.*T3.^2.*(( ...
  -1)+T.*((-1)+V)+(-1).*V).^(-1).*(1+T.*((-1)+V)+V)+2.*g.*r.*((1+g) ...
  .*T1+((-1)+g).*T2).*T3.*((-1)+T.*((-1)+V)+(-1).*V).^(-1).*(T.*(( ...
  -1)+V.^2)).^(1/2)+(-1).*r.*((1+g).*T1+((-1)+g).*T2).^2.*(1+T.*(( ...
  -1)+V)+V).*(4.*T.*((-1)+V)+(-4).*(1+V)).^(-1))).*pi.^(-1).*r.*(1+ ...
  T+V+(-1).*T.*V).^(-2).*((-1)+V.^2+(-2).*T.*((-3)+V.^2)+T.^2.*((-1) ...
  +V.^2)).^(-1).*((-1)+4.*T+10.*T.^2+4.*T.^3+(-1).*T.^4+(1/4).* ...
g.^2.*r.^4.*T.*((1+g).*T1+((-1)+g).*T2).^2.*T3.^2+(-1/2).*g.^2.* ...
 r.^4.*T.^2.*((1+g).*T1+((-1)+g).*T2).^2.*T3.^2+(1/4).*g.^2.*r.^4.* ...
  T.^3.*((1+g).*T1+((-1)+g).*T2).^2.*T3.^2+(-2).*V+12.*T.*V+(-12).* ...
  T.^3.*V+2.*T.^4.*V+8.*T.*V.^2+(-16).*T.^2.*V.^2+8.*T.^3.*V.^2+( ...
  -1/2).*g.^2.*r.^4.*T.*((1+g).*T1+((-1)+g).*T2).^2.*T3.^2.*V.^2+ ...
  g.^2.*r.^4.*T.^2.*((1+g).*T1+((-1)+g).*T2).^2.*T3.^2.*V.^2+(-1/2) ...
  .*g.^2.*r.^4.*T.^3.*((1+g).*T1+((-1)+g).*T2).^2.*T3.^2.*V.^2+2.* ...
  V.^3+(-4).*T.*V.^3+4.*T.^3.*V.^3+(-2).*T.^4.*V.^3+V.^4+(-4).*T.* ...
  V.^4+6.*T.^2.*V.^4+(-4).*T.^3.*V.^4+T.^4.*V.^4+(1/4).*g.^2.*r.^4.* ...
  T.*((1+g).*T1+((-1)+g).*T2).^2.*T3.^2.*V.^4+(-1/2).*g.^2.*r.^4.* ...
  T.^2.*((1+g).*T1+((-1)+g).*T2).^2.*T3.^2.*V.^4+(1/4).*g.^2.*r.^4.* ...
  T.^3.*((1+g).*T1+((-1)+g).*T2).^2.*T3.^2.*V.^4+g.^4.*r.^4.*((-1)+ ...
  T).^2.*T.*T3.^4.*((-1)+V.^2).^2+(-2).*g.*r.^2.*((1+g).*T1+((-1)+g) ...
  .*T2).*T3.*(T.*((-1)+V.^2)).^(1/2)+2.*g.*r.^2.*T.*((1+g).*T1+((-1) ...
  +g).*T2).*T3.*(T.*((-1)+V.^2)).^(1/2)+2.*g.*r.^2.*T.^2.*((1+g).* ...
  T1+((-1)+g).*T2).*T3.*(T.*((-1)+V.^2)).^(1/2)+(-2).*g.*r.^2.* ...
T.^3.*((1+g).*T1+((-1)+g).*T2).*T3.*(T.*((-1)+V.^2)).^(1/2)+(-3).* ...
  g.*r.^2.*((1+g).*T1+((-1)+g).*T2).*T3.*V.*(T.*((-1)+V.^2)).^(1/2)+ ...
  5.*g.*r.^2.*T.*((1+g).*T1+((-1)+g).*T2).*T3.*V.*(T.*((-1)+V.^2)) ...
  .^(1/2)+(-5).*g.*r.^2.*T.^2.*((1+g).*T1+((-1)+g).*T2).*T3.*V.*(T.* ...
  ((-1)+V.^2)).^(1/2)+3.*g.*r.^2.*T.^3.*((1+g).*T1+((-1)+g).*T2).* ...
  T3.*V.*(T.*((-1)+V.^2)).^(1/2)+g.*r.^2.*((1+g).*T1+((-1)+g).*T2).* ...
  T3.*V.^3.*(T.*((-1)+V.^2)).^(1/2)+(-3).*g.*r.^2.*T.*((1+g).*T1+(( ...
  -1)+g).*T2).*T3.*V.^3.*(T.*((-1)+V.^2)).^(1/2)+3.*g.*r.^2.*T.^2.*( ...
  (1+g).*T1+((-1)+g).*T2).*T3.*V.^3.*(T.*((-1)+V.^2)).^(1/2)+(-1).* ...
  g.*r.^2.*T.^3.*((1+g).*T1+((-1)+g).*T2).*T3.*V.^3.*(T.*((-1)+V.^2) ...
  ).^(1/2)+(1/4).*r.^4.*((-1)+T).^2.*((1+g).*T1+((-1)+g).*T2).^2.*(( ...
  -1)+V.^2).*((1/4).*T.*((1+g).*T1+((-1)+g).*T2).^2.*((-1)+V.^2)+ ...
  g.^2.*T.*T3.^2.*((-1)+V.^2)+(-1/2).*g.*((1+g).*T1+((-1)+g).*T2).* ...
  T3.*(1+T.*((-1)+V)+V).*(T.*((-1)+V.^2)).^(1/2))+g.^2.*r.^2.*((-1)+ ...
  T).*T3.^2.*((-1)+V.^2).*(T.^3.*((-1)+V).^2+(-1).*T.^2.*((-1)+V).*( ...
  (-5)+V+(1/2).*g.*r.^2.*((1+g).*T1+((-1)+g).*T2).*T3.*(T.*((-1)+ ...
  V.^2)).^(1/2))+(1+V).*(1+V+(1/2).*g.*r.^2.*((1+g).*T1+((-1)+g).* ...
  T2).*T3.*(T.*((-1)+V.^2)).^(1/2))+(-1).*T.*(5+6.*V+V.^2+g.*r.^2.*( ...
  (1+g).*T1+((-1)+g).*T2).*T3.*(T.*((-1)+V.^2)).^(1/2)))+(-1/2).* ...
  r.^2.*((-1)+T).*((1+g).*T1+((-1)+g).*T2).*((1/4).*g.*r.^2.*((-1)+ ...
  T).*T.^(-1).*((1+g).*T1+((-1)+g).*T2).^2.*T3.*(1+T.*((-1)+V)+V).*( ...
  T.*((-1)+V.^2)).^(3/2)+(-1/2).*((1+g).*T1+((-1)+g).*T2).*((-1)+ ...
  V.^2).*(T.^3.*((-1)+V).^2+(1+V).^2+g.^2.*r.^2.*((-1)+T).*T3.^2.*( ...
  1+T.*((-1)+V)+V).^2+(-1).*T.^2.*(5+(-6).*V+V.^2)+(-1).*T.*(5+6.*V+ ...
  V.^2))+g.*T3.*(T.*((-1)+V.^2)).^(1/2).*(g.^2.*r.^2.*((-1)+T).* ...
 T3.^2.*(1+T.*((-1)+V)+V).*((-1)+V.^2)+2.*(((-2)+V).*(1+V).^2+(-2) ...
  .*T.*V.*((-1)+V.^2)+T.^2.*(2+(-3).*V+V.^3))))))...
            ,0,inf,0,2*pi,-3*sigma,3*sigma,-3*sigma,3*sigma);
        
    case 'pspa'
        F=...
            integralN(@(r,theta,ra,ri)...
            exp(-(ra.^2 + ri.^2)./sigma).*1/(pi*sigma).*...
            (exp(1).^((-1/2).*r.*(r+2.*exp(1).^(sqrt(-1).*theta).*(ra+(sqrt(-1) ...
  *(-1)).*ri)+(-2).*exp(1).^((sqrt(-1)*(-1)).*theta).*(ra+sqrt(-1).* ...
  ri)+(-1).*exp(1).^(sqrt(-1).*theta).*(ra+(sqrt(-1)*(-1)).*ri).*( ...
  T1+g.*T1+T2+(-1).*g.*T2)+exp(1).^((sqrt(-1)*(-1)).*theta).*(ra+ ...
  sqrt(-1).*ri).*(T1+g.*T1+T2+(-1).*g.*T2)+(1/4).*r.*(T1+g.*T1+T2+( ...
  -1).*g.*T2).^2+(-1/2).*r.*((-2)+T1.^2+T2.^2+2.*g.*(T1.^2+(-1).* ...
  T2.^2)+g.^2.*((-4)+T1.^2+T2.^2+2.*T3.^2))+(-1).*g.^2.*r.*T3.^2.*(( ...
  -1)+T.*((-1)+V)+(-1).*V).^(-1).*(1+T.*((-1)+V)+V)+2.*g.*r.*((1+g) ...
  .*T1+((-1)+g).*T2).*T3.*((-1)+T.*((-1)+V)+(-1).*V).^(-1).*(T.*(( ...
  -1)+V.^2)).^(1/2)+(-1).*r.*((1+g).*T1+((-1)+g).*T2).^2.*(1+T.*(( ...
  -1)+V)+V).*(4.*T.*((-1)+V)+(-4).*(1+V)).^(-1))).*pi.^(-1).*r.*(1+ ...
  T.*((-1)+V)+V).^(-1).*(1+T+V+(-1).*T.*V).^(-2).*((1/4).*T.^3.*(4+ ...
  r.^2.*((1+g).*T1+((-1)+g).*T2).^2).*(1+g.^2.*r.^2.*T3.^2).*((-1)+ ...
  V).^3+(1+V).^2.*(1+V+g.*r.^2.*((1+g).*T1+((-1)+g).*T2).*T3.*(T.*(( ...
  -1)+V.^2)).^(1/2))+(-1).*T.*((-1)+V.^2).*(1+V+(-1).*g.*r.^2.*((1+ ...
  g).*T1+((-1)+g).*T2).*T3.*(T.*((-1)+V.^2)).^(1/2)+(1/8).*g.*r.^4.* ...
  ((1+g).*T1+((-1)+g).*T2).^3.*T3.*(T.*((-1)+V.^2)).^(1/2)+g.^2.* ...
  r.^2.*T3.^2.*(3+3.*V+(1/2).*g.*r.^2.*((1+g).*T1+((-1)+g).*T2).* ...
  T3.*(T.*((-1)+V.^2)).^(1/2))+(1/2).*r.^2.*((1+g).*T1+((-1)+g).*T2) ...
  .*((-1/2).*((1+g).*T1+((-1)+g).*T2).*((-3)+g.^2.*r.^2.*T3.^2).*(1+ ...
  V)+(1/4).*g.*r.^2.*((1+g).*T1+((-1)+g).*T2).^2.*T3.*(T.*((-1)+ ...
 V.^2)).^(1/2)+g.*T3.*((-2)+g.^2.*r.^2.*T3.^2).*(T.*((-1)+V.^2)).^( ...
  1/2)))+(1/16).*T.^2.*((-1)+V).^2.*((-16)+4.*g.^2.*r.^4.*((1+g).* ...
  T1+((-1)+g).*T2).^2.*T3.^2+(-16).*V+4.*g.^2.*r.^4.*((1+g).*T1+(( ...
  -1)+g).*T2).^2.*T3.^2.*V+16.*g.^4.*r.^4.*T3.^4.*(1+V)+(-24).*g.* ...
  r.^2.*((1+g).*T1+((-1)+g).*T2).*T3.*(T.*((-1)+V.^2)).^(1/2)+(-16) ...
  .*g.^2.*r.^2.*T3.^2.*((-2)+(-2).*V+(1/2).*g.*r.^2.*((1+g).*T1+(( ...
  -1)+g).*T2).*T3.*(T.*((-1)+V.^2)).^(1/2))+(-8).*r.^2.*((1+g).*T1+( ...
  (-1)+g).*T2).*((-1).*((1+g).*T1+((-1)+g).*T2).*(1+g.^2.*r.^2.* ...
 T3.^2).*(1+V)+(1/4).*g.*r.^2.*((1+g).*T1+((-1)+g).*T2).^2.*T3.*( ...
  T.*((-1)+V.^2)).^(1/2)+g.*T3.*(3+g.^2.*r.^2.*T3.^2).*(T.*((-1)+ ...
  V.^2)).^(1/2))+r.^4.*((1+g).*T1+((-1)+g).*T2).^2.*((1+g).^2.* ...
T1.^2.*(1+V)+((-1)+g).^2.*T2.^2.*(1+V)+4.*g.^2.*T3.^2.*(1+V)+(-2) ...
  .*((-1)+g).*g.*T2.*T3.*(T.*((-1)+V.^2)).^(1/2)+2.*(1+g).*T1.*((( ...
  -1)+g).*T2.*(1+V)+(-1).*g.*T3.*(T.*((-1)+V.^2)).^(1/2))))))...
            ,0,inf,0,2*pi,-3*sigma,3*sigma,-3*sigma,3*sigma);
        
    case 'paps'
        F=...
            integralN(@(r,theta,ra,ri)...
            exp(-(ra.^2 + ri.^2)./sigma).*1/(pi*sigma).*...
            ((1/4).*exp(1).^((-1/2).*r.*(r+2.*exp(1).^(sqrt(-1).*theta).*(ra+( ...
  sqrt(-1)*(-1)).*ri)+(-2).*exp(1).^((sqrt(-1)*(-1)).*theta).*(ra+ ...
  sqrt(-1).*ri)+(-1).*exp(1).^(sqrt(-1).*theta).*(ra+(sqrt(-1)*(-1)) ...
  .*ri).*(T1+g.*T1+T2+(-1).*g.*T2)+exp(1).^((sqrt(-1)*(-1)).*theta) ...
  .*(ra+sqrt(-1).*ri).*(T1+g.*T1+T2+(-1).*g.*T2)+(1/4).*r.*(T1+g.* ...
  T1+T2+(-1).*g.*T2).^2+(-1/2).*r.*((-2)+T1.^2+T2.^2+2.*g.*(T1.^2+( ...
  -1).*T2.^2)+g.^2.*((-4)+T1.^2+T2.^2+2.*T3.^2))+(-1).*g.^2.*r.* ...
 T3.^2.*((-1)+T.*((-1)+V)+(-1).*V).^(-1).*(1+T.*((-1)+V)+V)+2.*g.* ...
  r.*((1+g).*T1+((-1)+g).*T2).*T3.*((-1)+T.*((-1)+V)+(-1).*V).^(-1) ...
  .*(T.*((-1)+V.^2)).^(1/2)+(-1).*r.*((1+g).*T1+((-1)+g).*T2).^2.*( ...
  1+T.*((-1)+V)+V).*(4.*T.*((-1)+V)+(-4).*(1+V)).^(-1))).*pi.^(-1).* ...
  r.*(1+T.*((-1)+V)+V).^(-1).*(1+T+V+(-1).*T.*V).^(-2).*(4+4.*T+(-4) ...
  .*T.^2+(-4).*T.^3+(-1).*g.^2.*r.^4.*T.*((1+g).*T1+((-1)+g).*T2) ...
  .^2.*T3.^2+12.*V+4.*T.*V+4.*T.^2.*V+12.*T.^3.*V+(-1).*g.^2.*r.^4.* ...
  T.*((1+g).*T1+((-1)+g).*T2).^2.*T3.^2.*V+12.*V.^2+(-4).*T.*V.^2+ ...
  4.*T.^2.*V.^2+(-12).*T.^3.*V.^2+g.^2.*r.^4.*T.*((1+g).*T1+((-1)+g) ...
  .*T2).^2.*T3.^2.*V.^2+4.*V.^3+(-4).*T.*V.^3+(-4).*T.^2.*V.^3+4.* ...
  T.^3.*V.^3+g.^2.*r.^4.*T.*((1+g).*T1+((-1)+g).*T2).^2.*T3.^2.* ...
 V.^3+4.*g.^4.*r.^4.*T.*T3.^4.*((-1)+V).*(1+V).^2+6.*g.*r.^2.*((1+ ...
  g).*T1+((-1)+g).*T2).*T3.*(T.*((-1)+V.^2)).^(1/2)+4.*g.*r.^2.*T.*( ...
  (1+g).*T1+((-1)+g).*T2).*T3.*(T.*((-1)+V.^2)).^(1/2)+(-2).*g.* ...
 r.^2.*T.^2.*((1+g).*T1+((-1)+g).*T2).*T3.*(T.*((-1)+V.^2)).^(1/2)+ ...
  12.*g.*r.^2.*((1+g).*T1+((-1)+g).*T2).*T3.*V.*(T.*((-1)+V.^2)).^( ...
  1/2)+4.*g.*r.^2.*T.^2.*((1+g).*T1+((-1)+g).*T2).*T3.*V.*(T.*((-1)+ ...
  V.^2)).^(1/2)+6.*g.*r.^2.*((1+g).*T1+((-1)+g).*T2).*T3.*V.^2.*(T.* ...
  ((-1)+V.^2)).^(1/2)+(-4).*g.*r.^2.*T.*((1+g).*T1+((-1)+g).*T2).* ...
  T3.*V.^2.*(T.*((-1)+V.^2)).^(1/2)+(-2).*g.*r.^2.*T.^2.*((1+g).*T1+ ...
  ((-1)+g).*T2).*T3.*V.^2.*(T.*((-1)+V.^2)).^(1/2)+r.^4.*((1+g).*T1+ ...
  ((-1)+g).*T2).^2.*(1+V).*((1/4).*T.*((1+g).*T1+((-1)+g).*T2).^2.*( ...
  (-1)+V.^2)+g.^2.*T.*T3.^2.*((-1)+V.^2)+(-1/2).*g.*((1+g).*T1+((-1) ...
  +g).*T2).*T3.*(1+T.*((-1)+V)+V).*(T.*((-1)+V.^2)).^(1/2))+(-2).* ...
  r.^2.*((1+g).*T1+((-1)+g).*T2).*((1/4).*g.*r.^2.*((1+g).*T1+((-1)+ ...
  g).*T2).^2.*T3.*(1+V).*(1+T.*((-1)+V)+V).*(T.*((-1)+V.^2)).^(1/2)+ ...
  (-1/2).*((1+g).*T1+((-1)+g).*T2).*(1+V).*(3.*T.^2.*((-1)+V).^2+( ...
  -1).*(1+V).^2+g.^2.*r.^2.*T3.^2.*(1+T.*((-1)+V)+V).^2+(-2).*T.*(( ...
  -1)+V.^2))+g.*T3.*(T.*((-1)+V.^2)).^(1/2).*(T.^2.*((-1)+V).^2+(-3) ...
  .*(1+V).^2+g.^2.*r.^2.*T3.^2.*(1+V).*(1+T.*((-1)+V)+V)+2.*T.*((-1) ...
  +V.^2)))+(-4).*g.^2.*r.^2.*T3.^2.*(1+V).*((-3).*T.^2.*((-1)+V).^2+ ...
  (1+V).*(1+V+(1/2).*g.*r.^2.*((1+g).*T1+((-1)+g).*T2).*T3.*(T.*(( ...
  -1)+V.^2)).^(1/2))+T.*((-1)+V).*(2+2.*V+(1/2).*g.*r.^2.*((1+g).* ...
  T1+((-1)+g).*T2).*T3.*(T.*((-1)+V.^2)).^(1/2)))))...
            ,0,inf,0,2*pi,-3*sigma,3*sigma,-3*sigma,3*sigma);
end

F = real(F);

end