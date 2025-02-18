function [F] = sq_loss_eta(T, l, t, e, operation, g, eta, ra)
warning('off')

switch operation
    case 'epr'
        % send mode B
        F=...
            integral2(@(r,theta)...
            (exp(1).^((-1/4).*r.^2.*((-2).*((-1)+l.^2).^(-1).*(1+l.^2+(-4).* ...
  eta.*g.*l.*t.^(1/2)+g.^2.*(1+e.*eta.^2.*((-1)+l.^2).*((-1)+t)+ ...
  l.^2.*((-1)+2.*eta.^2.*t)))+2.*(1+eta.^2.*g.^2).*cosh(2.*ra)+exp( ...
  1).^((sqrt(-1)*(-2)).*theta).*(1+exp(1).^((sqrt(-1)*4).*theta)).*( ...
  1+eta.^2.*g.^2).*sinh(2.*ra))).*pi.^(-1).*r)...
            ,0,inf,0,2*pi);
    case 'ps'
        F=...
            integral2(@(r,theta)...
            (exp(1).^((-1/4).*r.^2.*((-2).*((-1)+eta.^2).*g.^2+(-2).*((-2)+ ...
l.^2.*((-2)+((-2)+e).*t.*((-1)+T))+e.*t.*((-1)+T)).*(2+l.^2.*((-2) ...
  +((-2)+e).*t.*((-1)+T))+e.*(t+(-1).*t.*T)).^(-1)+(-16).*eta.*g.* ...
  l.*t.^(1/2).*T.^(1/2).*(2+l.^2.*((-2)+((-2)+e).*t.*((-1)+T))+e.*( ...
  t+(-1).*t.*T)).^(-1)+(-2).*eta.^2.*g.^2.*(2+l.^2.*((-2)+((-2)+e).* ...
  t.*((-1)+T))+e.*(t+(-1).*t.*T)).^(-1).*((-2)+(-1).*e.*t.*(1+T)+ ...
  l.^2.*(2+((-2)+e).*t.*(1+T)))+2.*cosh(2.*ra)+2.*eta.^2.*g.^2.* ...
 cosh(2.*ra)+exp(1).^((sqrt(-1)*(-2)).*theta).*(1+exp(1).^((sqrt( ...
  -1)*4).*theta)).*sinh(2.*ra)+exp(1).^((sqrt(-1)*(-2)).*theta).*(1+ ...
  exp(1).^((sqrt(-1)*4).*theta)).*eta.^2.*g.^2.*sinh(2.*ra))).*((-2) ...
  .*l.^2+e.*((-1)+l.^2)).^(-1).*pi.^(-1).*r.*(2+l.^2.*((-2)+((-2)+e) ...
  .*t.*((-1)+T))+e.*(t+(-1).*t.*T)).^(-1).*(e.^2.*((-1)+l.^2).^2.* ...
  t.*((-1)+T+eta.^2.*g.^2.*r.^2.*T)+4.*l.^2.*((-1)+r.^2+(-2).*eta.* ...
  g.*l.*r.^2.*t.^(1/2).*T.^(1/2)+l.^2.*(1+t.*((-1)+T+eta.^2.*g.^2.* ...
  r.^2.*T)))+(-2).*e.*((-1)+l.^2).*((-1)+(-2).*eta.*g.*l.*r.^2.*t.^( ...
  1/2).*T.^(1/2)+l.^2.*(1+2.*t.*((-1)+T+eta.^2.*g.^2.*r.^2.*T)))))...
            ,0,inf,0,2*pi);
    case 'pa'
        F=...
            integral2(@(r,theta)....
            (exp(1).^((-1/4).*r.^2.*((-2).*((-1)+eta.^2).*g.^2+(-2).*((-2)+ ...
 l.^2.*((-2)+((-2)+e).*t.*((-1)+T))+e.*t.*((-1)+T)).*(2+l.^2.*((-2) ...
  +((-2)+e).*t.*((-1)+T))+e.*(t+(-1).*t.*T)).^(-1)+(-16).*eta.*g.* ...
  l.*t.^(1/2).*T.^(1/2).*(2+l.^2.*((-2)+((-2)+e).*t.*((-1)+T))+e.*( ...
  t+(-1).*t.*T)).^(-1)+(-2).*eta.^2.*g.^2.*(2+l.^2.*((-2)+((-2)+e).* ...
  t.*((-1)+T))+e.*(t+(-1).*t.*T)).^(-1).*((-2)+(-1).*e.*t.*(1+T)+ ...
  l.^2.*(2+((-2)+e).*t.*(1+T)))+2.*cosh(2.*ra)+2.*eta.^2.*g.^2.* ...
 cosh(2.*ra)+exp(1).^((sqrt(-1)*(-2)).*theta).*(1+exp(1).^((sqrt( ...
  -1)*4).*theta)).*sinh(2.*ra)+exp(1).^((sqrt(-1)*(-2)).*theta).*(1+ ...
  exp(1).^((sqrt(-1)*4).*theta)).*eta.^2.*g.^2.*sinh(2.*ra))).*pi.^( ...
  -1).*r.*((-2)+(-1).*e.*t+l.^2.*(2+((-2)+e).*t)).^(-1).*(2+l.^2.*(( ...
  -2)+((-2)+e).*t.*((-1)+T))+e.*(t+(-1).*t.*T)).^(-1).*(l.^4.*(2+(( ...
  -2)+e).*t).*((-2)+eta.^2.*g.^2.*r.^2.*(2+((-2)+e).*t)+((-2)+e).* ...
  t.*((-1)+T))+(2+e.*t).*((-2)+eta.^2.*g.^2.*r.^2.*(2+e.*t)+e.*t.*(( ...
  -1)+T))+4.*eta.*g.*l.^3.*r.^2.*t.^(1/2).*(2+((-2)+e).*t).*T.^(1/2) ...
  +(-4).*eta.*g.*l.*r.^2.*t.^(1/2).*(2+e.*t).*T.^(1/2)+l.^2.*(8+(-2) ...
  .*eta.^2.*g.^2.*r.^2.*(4+4.*((-1)+e).*t+((-2)+e).*e.*t.^2)+(-2).*( ...
  (-2)+e).*e.*t.^2.*((-1)+T)+4.*t.*((-2)+2.*e+T+(-1).*e.*T+r.^2.*T)) ...
  ))...
            ,0,inf,0,2*pi);
    case 'pc'
        F=...
            integral2(@(r,theta)...
            (exp(1).^((-1/4).*r.^2.*((-2).*((-1)+eta.^2).*g.^2+(-2).*((-2)+ ...
l.^2.*((-2)+((-2)+e).*t.*((-1)+T))+e.*t.*((-1)+T)).*(2+l.^2.*((-2) ...
  +((-2)+e).*t.*((-1)+T))+e.*(t+(-1).*t.*T)).^(-1)+(-16).*eta.*g.* ...
  l.*t.^(1/2).*T.^(1/2).*(2+l.^2.*((-2)+((-2)+e).*t.*((-1)+T))+e.*( ...
  t+(-1).*t.*T)).^(-1)+(-2).*eta.^2.*g.^2.*(2+l.^2.*((-2)+((-2)+e).* ...
  t.*((-1)+T))+e.*(t+(-1).*t.*T)).^(-1).*((-2)+(-1).*e.*t.*(1+T)+ ...
  l.^2.*(2+((-2)+e).*t.*(1+T)))+2.*cosh(2.*ra)+2.*eta.^2.*g.^2.* ...
 cosh(2.*ra)+exp(1).^((sqrt(-1)*(-2)).*theta).*(1+exp(1).^((sqrt( ...
  -1)*4).*theta)).*sinh(2.*ra)+exp(1).^((sqrt(-1)*(-2)).*theta).*(1+ ...
  exp(1).^((sqrt(-1)*4).*theta)).*eta.^2.*g.^2.*sinh(2.*ra))).*pi.^( ...
  -1).*r.*(2+l.^2.*((-2)+((-2)+e).*t.*((-1)+T))+e.*(t+(-1).*t.*T)) ...
  .^(-2).*(e.^2.*((-1)+l.^2).^2.*t.^2.*((-1)+T).^2+(-2).*e.*((-1)+ ...
  l.^2).*t.*(1+l.^2.*((-1)+2.*t)).*((-1)+T).^2+4.*(l.^2.*(t.*((-1)+ ...
  T).^2+(-2).*T)+T+l.^4.*((-1).*t.*((-1)+T).^2+t.^2.*((-1)+T).^2+T)) ...
  ).^(-1).*(e.^4.*((-1)+l.^2).^4.*t.^4.*((-1)+T).^2.*(((-1)+T).^2+ ...
  eta.^4.*g.^4.*r.^4.*T+eta.^2.*g.^2.*r.^2.*((-1)+T.^2))+(-2).* ...
e.^3.*((-1)+l.^2).^3.*t.^3.*((-1)+T).^2.*(3+(-4).*T+2.*eta.^4.* ...
g.^4.*r.^4.*T+T.^2+eta.^2.*g.^2.*r.^2.*((-3)+2.*T+T.^2)+(-2).* ...
eta.*g.*l.*r.^2.*t.^(1/2).*T.^(1/2).*(2.*((-1)+T)+eta.^2.*g.^2.* ...
 r.^2.*(1+T))+l.^2.*((3+4.*t.*((-1)+T)+(-1).*T).*((-1)+T)+2.* ...
eta.^4.*g.^4.*r.^4.*((-1)+2.*t).*T+eta.^2.*g.^2.*r.^2.*((-1)+T).*( ...
  (-3)+(-1).*T+4.*t.*(1+T))))+4.*e.^2.*((-1)+l.^2).^2.*t.^2.*((-1)+ ...
  T).*((-3)+3.*eta.^2.*g.^2.*r.^2+4.*T+(-7).*eta.^2.*g.^2.*r.^2.*T+( ...
  -1).*eta.^4.*g.^4.*r.^4.*T+(-1).*T.^2+2.*eta.^2.*g.^2.*r.^2.*T.^2+ ...
  eta.^4.*g.^4.*r.^4.*T.^2+(-2).*eta.*g.*l.*r.^2.*t.^(1/2).*((-1)+T) ...
  .*T.^(1/2).*(3.*((-1)+T)+eta.^2.*g.^2.*r.^2.*(2+T))+(-2).*eta.*g.* ...
  l.^3.*r.^2.*t.^(1/2).*((-1)+T).*T.^(1/2).*(3.*((-1)+2.*t).*((-1)+ ...
  T)+eta.^2.*g.^2.*r.^2.*((-2)+(-1).*T+3.*t.*(1+T)))+l.^4.*(eta.^4.* ...
  g.^4.*r.^4.*(1+(-6).*t+6.*t.^2).*((-1)+T).*T+eta.^2.*g.^2.*r.^2.*( ...
  3+(-7).*T+2.*T.^2+6.*t.^2.*((-1)+T).^2.*(1+T)+(-3).*t.*((-1)+T) ...
  .^2.*(3+T))+((-1)+T).*(3+6.*t.^2.*((-1)+T).^2+(-1).*T+(-3).*t.*(3+ ...
  (-4).*T+T.^2)))+l.^2.*(2.*eta.^4.*g.^4.*r.^4.*((-1)+3.*t).*((-1)+ ...
  T).*T+((-1)+T).*(2.*((-3)+T)+t.*((-1)+T).*(3.*((-3)+T)+r.^2.*(1+T) ...
  ))+eta.^2.*g.^2.*r.^2.*((-6)+14.*T+(-4).*T.^2+t.*((-1)+T).*(3.*(( ...
  -3)+2.*T+T.^2)+r.^2.*(1+4.*T+T.^2)))))+16.*(2.*eta.*g.*l.*r.^2.* ...
  t.^(1/2).*((-1)+T).*T.^(1/2)+T+(-2).*eta.*g.*l.^3.*r.^2.*t.^(1/2) ...
  .*((-1)+T).*T.^(1/2).*(3+(1+eta.^2.*g.^2).*r.^2.*t.*((-1)+T)+(-2) ...
  .*t.*T)+(-2).*eta.*g.*l.^7.*r.^2.*t.^(1/2).*((-1)+T).*T.^(1/2).*( ...
  1+(-3).*t.^2+2.*t.^3+(-2).*t.*T+6.*t.^2.*T+(-4).*t.^3.*T+(-3).* ...
  t.^2.*T.^2+2.*t.^3.*T.^2+eta.^2.*g.^2.*r.^2.*((-1)+t).*t.*((-1)+T) ...
  .*((-1)+t+t.*T))+l.^2.*((-4).*T+t.*((-1)+T).*((-1)+eta.^2.*g.^2.* ...
  r.^4.*((-1)+T)+(-1).*T+(-1).*(1+eta.^2.*g.^2).*r.^2.*((-1)+3.*T))) ...
  +l.^8.*(T+t.*((-1)+T).*(1+T+eta.^2.*g.^2.*r.^2.*((-1)+3.*T))+ ...
t.^4.*((-1)+T).^2.*(((-1)+T).^2+eta.^4.*g.^4.*r.^4.*T+eta.^2.* ...
g.^2.*r.^2.*((-1)+T.^2))+(-1).*t.^3.*((-1)+T).^2.*(3+(-4).*T+2.* ...
 eta.^4.*g.^4.*r.^4.*T+T.^2+eta.^2.*g.^2.*r.^2.*((-3)+2.*T+T.^2))+ ...
  t.^2.*((-1).*((-3)+T).*((-1)+T).^2+eta.^4.*g.^4.*r.^4.*((-1)+T) ...
  .^2.*T+eta.^2.*g.^2.*r.^2.*((-3)+10.*T+(-9).*T.^2+2.*T.^3)))+(-2) ...
  .*eta.*g.*l.^5.*r.^2.*t.^(1/2).*((-1)+T).*T.^(1/2).*((-3)+3.*t.^2+ ...
  4.*t.*T+(-6).*t.^2.*T+3.*t.^2.*T.^2+r.^2.*t.*((-1)+T).*((-1)+t+t.* ...
  T+eta.^2.*g.^2.*((-2)+t.*(2+T))))+l.^6.*((-4).*T+t.*((-1)+T).*( ...
  eta.^2.*g.^2.*r.^4.*((-1)+T)+(-3).*(1+T)+(-1).*(1+3.*eta.^2.*g.^2) ...
  .*r.^2.*((-1)+3.*T))+t.^3.*((-1)+T).^2.*(3+(-4).*T+T.^2+eta.^2.* ...
  g.^2.*r.^4.*(1+2.*(2+eta.^2.*g.^2).*T+T.^2)+r.^2.*((-1)+T).*(1+T+ ...
  eta.^2.*g.^2.*(3+T)))+(-2).*t.^2.*((-1)+T).*((-3)+4.*T+(-1).*T.^2+ ...
  eta.^2.*g.^2.*r.^4.*((-1)+T).*(1+(2+eta.^2.*g.^2).*T)+r.^2.*(((-1) ...
  +T).^2+eta.^2.*g.^2.*(3+(-7).*T+2.*T.^2))))+l.^4.*(6.*T+(-1).*t.*( ...
  (-1)+T).*(2.*eta.^2.*g.^2.*r.^4.*((-1)+T)+(-3).*(1+T)+(-1).*(2+3.* ...
  eta.^2.*g.^2).*r.^2.*((-1)+3.*T))+t.^2.*((-1)+T).*((-1).*((-3)+T) ...
  .*((-1)+T)+r.^4.*((-1)+T).*(T+eta.^4.*g.^4.*T+2.*eta.^2.*g.^2.*(1+ ...
  2.*T))+r.^2.*(2.*((-1)+T).^2+eta.^2.*g.^2.*(3+(-7).*T+2.*T.^2))))) ...
  +(-8).*e.*((-1)+l.^2).*t.*((-1)+T).*((-1)+eta.^2.*g.^2.*r.^2+(-1) ...
  .*T+(-3).*eta.^2.*g.^2.*r.^2.*T+2.*eta.*g.*l.*r.^2.*t.^(1/2).*T.^( ...
  1/2).*((-1).*eta.^2.*g.^2.*r.^2.*((-1)+T)+2.*T)+(-2).*eta.*g.* ...
 l.^5.*r.^2.*t.^(1/2).*T.^(1/2).*(2.*((-3).*t.*((-1)+T).^2+3.* ...
t.^2.*((-1)+T).^2+(-1).*T)+eta.^2.*g.^2.*r.^2.*((-1)+T).*(1+3.* ...
t.^2.*(1+T)+(-2).*t.*(2+T)))+(-2).*eta.*g.*l.^3.*r.^2.*t.^(1/2).* ...
  T.^(1/2).*(4.*T+t.*((-1)+T).*(6.*((-1)+T)+r.^2.*(1+T))+2.*eta.^2.* ...
  g.^2.*r.^2.*((-1)+T).*((-1)+t.*(2+T)))+l.^6.*(1+(-6).*t+9.*t.^2+( ...
  -4).*t.^3+T+8.*t.*T+(-21).*t.^2.*T+12.*t.^3.*T+2.*eta.^4.*g.^4.* ...
  r.^4.*t.*(1+(-3).*t+2.*t.^2).*((-1)+T).*T+(-2).*t.*T.^2+15.*t.^2.* ...
  T.^2+(-12).*t.^3.*T.^2+(-3).*t.^2.*T.^3+4.*t.^3.*T.^3+eta.^2.* ...
 g.^2.*r.^2.*((-1)+3.*T+4.*t.^3.*((-1)+T).^2.*(1+T)+(-3).*t.^2.*(( ...
 -1)+T).^2.*(3+T)+2.*t.*(3+(-7).*T+2.*T.^2)))+l.^4.*((-2).*t.*(6+ ...
  r.^2.*((-1)+T)+(-2).*T).*((-1)+T)+2.*eta.^4.*g.^4.*r.^4.*t.*((-2)+ ...
  3.*t).*((-1)+T).*T+(-3).*(1+T)+t.^2.*((-1)+T).^2.*(3.*((-3)+T)+2.* ...
  r.^2.*(1+T))+eta.^2.*g.^2.*r.^2.*(3+(-9).*T+2.*t.*((-6)+14.*T+(-4) ...
  .*T.^2+r.^2.*(1+T+(-2).*T.^2))+t.^2.*((-1)+T).*(3.*((-3)+2.*T+ ...
 T.^2)+2.*r.^2.*(1+4.*T+T.^2))))+l.^2.*(2.*t.*(3+r.^2.*((-1)+T)+( ...
  -1).*T).*((-1)+T)+2.*eta.^4.*g.^4.*r.^4.*t.*((-1)+T).*T+3.*(1+T)+ ...
  eta.^2.*g.^2.*r.^2.*((-3)+9.*T+2.*t.*(3+(-7).*T+2.*T.^2+r.^2.*(( ...
  -1)+(-1).*T+2.*T.^2)))))))...
            ,0,inf,0,2*pi);
        
    case 'sa'
        F=...
            integral2(@(r,theta)...
            (exp(1).^((1/2).*(((-1)+eta.^2).*g.^2.*r.^2+8.*eta.*g.*l.*r.^2.* ...
  t.^(1/2).*T.*(2+e.*(t+(-1).*t.*T.^2)+l.^2.*((-2)+((-2)+e).*t.*(( ...
  -1)+T.^2))).^(-1)+r.^2.*((-2)+e.*t.*((-1)+T.^2)+l.^2.*((-2)+((-2)+ ...
  e).*t.*((-1)+T.^2))).*(2+e.*(t+(-1).*t.*T.^2)+l.^2.*((-2)+((-2)+e) ...
  .*t.*((-1)+T.^2))).^(-1)+eta.^2.*g.^2.*r.^2.*(2+e.*(t+(-1).*t.* ...
  T.^2)+l.^2.*((-2)+((-2)+e).*t.*((-1)+T.^2))).^(-1).*((-2)+(-1).* ...
  e.*t.*(1+T.^2)+l.^2.*(2+((-2)+e).*t.*(1+T.^2)))+(-1).*abs(exp(1) ...
  .^((sqrt(-1)*(-1)).*theta).*r.*(exp(1).^((sqrt(-1)*2).*theta).* ...
  cosh(ra)+sinh(ra))).^2+(-1).*abs((-1).*exp(1).^((sqrt(-1)*(-1)).* ...
  theta).*eta.*g.*r.*(exp(1).^((sqrt(-1)*2).*theta).*cosh(ra)+sinh( ...
  ra))).^2)).*((-2).*l.^2+e.*((-1)+l.^2)).^(-1).*pi.^(-1).*r.*(2+e.* ...
  (t+(-1).*t.*T.^2)+l.^2.*((-2)+((-2)+e).*t.*((-1)+T.^2))).^(-2).*(( ...
  -2)+(-1).*e.*t.*(1+T.^2)+l.^2.*(2+((-2)+e).*t.*(1+T.^2))).^(-1).*( ...
  e.^4.*((-1)+l.^2).^4.*t.^3.*(eta.^4.*g.^4.*r.^4.*T.^2+((-1)+T.^2) ...
  .^2.*(1+T.^2)+eta.^2.*g.^2.*r.^2.*((-1)+(-2).*T.^2+3.*T.^4))+(-2) ...
  .*e.^3.*((-1)+l.^2).^3.*t.^2.*(3+(-2).*T.^2+2.*eta.^4.*g.^4.* ...
  r.^4.*T.^2+(-1).*T.^4+eta.^2.*g.^2.*r.^2.*((-3)+(-4).*T.^2+3.* ...
  T.^4)+(-2).*eta.*g.*l.*r.^2.*t.^(1/2).*T.*((-3)+2.*T.^2+T.^4+ ...
  eta.^2.*g.^2.*r.^2.*(1+T.^2))+l.^2.*(2.*eta.^4.*g.^4.*r.^4.*((-1)+ ...
  2.*t).*T.^2+((-1)+T.^2).*(3+T.^2+4.*t.*((-1)+T.^4))+eta.^2.*g.^2.* ...
  r.^2.*(3+4.*T.^2+(-3).*T.^4+4.*t.*((-1)+(-2).*T.^2+3.*T.^4))))+4.* ...
  e.^2.*((-1)+l.^2).^2.*t.*(3+(-3).*eta.^2.*g.^2.*r.^2+(-1).*T.^2+( ...
  -2).*eta.^2.*g.^2.*r.^2.*T.^2+eta.^4.*g.^4.*r.^4.*T.^2+(-2).*eta.* ...
  g.*l.*r.^2.*t.^(1/2).*T.*(2.*((-3)+T.^2)+eta.^2.*g.^2.*r.^2.*(2+ ...
  T.^2))+(-2).*eta.*g.*l.^3.*r.^2.*t.^(1/2).*T.*(6+(-2).*T.^2+3.*t.* ...
  ((-3)+2.*T.^2+T.^4)+eta.^2.*g.^2.*r.^2.*((-2)+(-1).*T.^2+3.*t.*(1+ ...
  T.^2)))+l.^4.*(3+(-1).*T.^2+eta.^4.*g.^4.*r.^4.*(1+(-6).*t+6.* ...
  t.^2).*T.^2+6.*t.^2.*((-1)+T.^2).^2.*(1+T.^2)+3.*t.*((-3)+2.*T.^2+ ...
  T.^4)+eta.^2.*g.^2.*r.^2.*((-3)+(-2).*T.^2+t.*(9+12.*T.^2+(-9).* ...
  T.^4)+6.*t.^2.*((-1)+(-2).*T.^2+3.*T.^4)))+l.^2.*(2.*eta.^4.* ...
  g.^4.*r.^4.*((-1)+3.*t).*T.^2+2.*((-3)+T.^2)+t.*((-1)+T.^2).*((-9) ...
  +r.^2+(-3).*T.^2+3.*r.^2.*T.^2)+eta.^2.*g.^2.*r.^2.*(6+4.*T.^2+t.* ...
  ((-9)+(-12).*T.^2+9.*T.^4+r.^2.*(1+4.*T.^2+T.^4)))))+(-8).*e.*(( ...
  -1)+l.^2).*(1+(-1).*eta.^2.*g.^2.*r.^2+(-2).*eta.*g.*l.*r.^2.*(( ...
  -3)+eta.^2.*g.^2.*r.^2).*t.^(1/2).*T+(-2).*eta.*g.*l.^5.*r.^2.* ...
  t.^(1/2).*T.*((-3)+(-4).*t.*((-3)+T.^2)+3.*t.^2.*((-3)+2.*T.^2+ ...
  T.^4)+eta.^2.*g.^2.*r.^2.*(1+3.*t.^2.*(1+T.^2)+(-2).*t.*(2+T.^2))) ...
  +(-2).*eta.*g.*l.^3.*r.^2.*t.^(1/2).*T.*(6+4.*t.*((-3)+T.^2)+ ...
  r.^2.*t.*(1+T.^2)+2.*eta.^2.*g.^2.*r.^2.*((-1)+t.*(2+T.^2)))+ ...
  l.^2.*((-3)+2.*eta.^4.*g.^4.*r.^4.*t.*T.^2+(-2).*t.*((-3)+T.^2+ ...
  r.^2.*(1+T.^2))+eta.^2.*g.^2.*r.^2.*(3+2.*t.*((-3)+r.^2+(-2).* ...
  T.^2+2.*r.^2.*T.^2)))+l.^6.*((-1)+2.*eta.^4.*g.^4.*r.^4.*t.*(1+( ...
  -3).*t+2.*t.^2).*T.^2+(-2).*t.*((-3)+T.^2)+4.*t.^3.*((-1)+T.^2) ...
  .^2.*(1+T.^2)+3.*t.^2.*((-3)+2.*T.^2+T.^4)+eta.^2.*g.^2.*r.^2.*(1+ ...
  (-2).*t.*(3+2.*T.^2)+t.^2.*(9+12.*T.^2+(-9).*T.^4)+4.*t.^3.*((-1)+ ...
  (-2).*T.^2+3.*T.^4)))+l.^4.*(3+2.*eta.^4.*g.^4.*r.^4.*t.*((-2)+3.* ...
  t).*T.^2+2.*t.*(2.*((-3)+T.^2)+r.^2.*(1+T.^2))+t.^2.*((-1)+T.^2).* ...
  ((-3).*(3+T.^2)+r.^2.*(2+6.*T.^2))+eta.^2.*g.^2.*r.^2.*((-3)+(-2) ...
  .*t.*((-6)+(-4).*T.^2+r.^2.*(1+2.*T.^2))+t.^2.*((-9)+(-12).*T.^2+ ...
  9.*T.^4+2.*r.^2.*(1+4.*T.^2+T.^4)))))+16.*l.^2.*(((-1)+r.^2).*(( ...
  -1)+eta.^2.*g.^2.*r.^2)+(-2).*eta.*g.*l.*r.^2.*((-3)+(1+eta.^2.* ...
  g.^2).*r.^2).*t.^(1/2).*T+(-2).*eta.*g.*l.^5.*r.^2.*t.^(1/2).*T.*( ...
  (-3)+(-2).*t.*((-3)+T.^2)+t.^2.*((-3)+2.*T.^2+T.^4)+eta.^2.*g.^2.* ...
  r.^2.*(1+t.^2.*(1+T.^2)+(-1).*t.*(2+T.^2)))+l.^2.*((-3)+(-1).*t.*( ...
  (-3)+T.^2)+r.^4.*(t.*T.^2+eta.^4.*g.^4.*t.*T.^2+2.*eta.^2.*g.^2.*( ...
  (-1)+t+2.*t.*T.^2))+(-1).*r.^2.*(2.*((-1)+t+t.*T.^2)+eta.^2.* ...
  g.^2.*((-3)+3.*t+2.*t.*T.^2)))+l.^6.*(eta.^4.*g.^4.*r.^4.*((-1)+t) ...
  .^2.*t.*T.^2+((-1)+t+t.*T.^2).*(1+t.*((-1)+T.^2)).^2+eta.^2.* ...
  g.^2.*r.^2.*(1+(-1).*t.*(3+2.*T.^2)+t.^2.*(3+4.*T.^2+(-3).*T.^4)+ ...
  t.^3.*((-1)+(-2).*T.^2+3.*T.^4)))+(-2).*eta.*g.*l.^3.*r.^2.*t.^( ...
  1/2).*T.*(6+2.*t.*((-3)+T.^2)+r.^2.*((-1)+t+t.*T.^2+eta.^2.*g.^2.* ...
  ((-2)+t.*(2+T.^2))))+l.^4.*(3+2.*t.*((-3)+T.^2)+(-1).*t.^2.*((-3)+ ...
  2.*T.^2+T.^4)+eta.^2.*g.^2.*r.^4.*(1+(-2).*t.*(1+(2+eta.^2.*g.^2) ...
  .*T.^2)+t.^2.*(1+2.*(2+eta.^2.*g.^2).*T.^2+T.^4))+r.^2.*((-1)+2.* ...
  t.*(1+T.^2)+t.^2.*((-1)+(-2).*T.^2+3.*T.^4)+eta.^2.*g.^2.*((-3)+ ...
  t.*(6+4.*T.^2)+t.^2.*((-3)+(-4).*T.^2+3.*T.^4)))))))...
            ,0,inf,0,2*pi);
        
    case 'as'
        F=...
            integral2(@(r,theta)...
            (exp(1).^((1/2).*(((-1)+eta.^2).*g.^2.*r.^2+8.*eta.*g.*l.*r.^2.* ...
t.^(1/2).*T.*(2+e.*(t+(-1).*t.*T.^2)+l.^2.*((-2)+((-2)+e).*t.*(( ...
 -1)+T.^2))).^(-1)+r.^2.*((-2)+e.*t.*((-1)+T.^2)+l.^2.*((-2)+((-2)+ ...
  e).*t.*((-1)+T.^2))).*(2+e.*(t+(-1).*t.*T.^2)+l.^2.*((-2)+((-2)+e) ...
  .*t.*((-1)+T.^2))).^(-1)+eta.^2.*g.^2.*r.^2.*(2+e.*(t+(-1).*t.* ...
  T.^2)+l.^2.*((-2)+((-2)+e).*t.*((-1)+T.^2))).^(-1).*((-2)+(-1).* ...
  e.*t.*(1+T.^2)+l.^2.*(2+((-2)+e).*t.*(1+T.^2)))+(-1).*abs(exp(1) ...
  .^((sqrt(-1)*(-1)).*theta).*r.*(exp(1).^((sqrt(-1)*2).*theta).* ...
  cosh(ra)+sinh(ra))).^2+(-1).*abs((-1).*exp(1).^((sqrt(-1)*(-1)).* ...
  theta).*eta.*g.*r.*(exp(1).^((sqrt(-1)*2).*theta).*cosh(ra)+sinh( ...
  ra))).^2)).*pi.^(-1).*r.*((-2)+(-1).*e.*t+l.^2.*(2+((-2)+e).*t)) ...
  .^(-1).*(2+e.*(t+(-1).*t.*T.^2)+l.^2.*((-2)+((-2)+e).*t.*((-1)+ ...
  T.^2))).^(-2).*((-2)+(-1).*e.*t.*(1+T.^2)+l.^2.*(2+((-2)+e).*t.*( ...
  1+T.^2))).^(-1).*(16+(-8).*e.*t.*((-4)+(1+3.*eta.^2.*g.^2.*r.^2).* ...
  T.^2)+4.*e.^2.*t.^2.*(6+(-3).*(1+3.*eta.^2.*g.^2.*r.^2).*T.^2+(( ...
  -1)+2.*eta.^2.*g.^2.*r.^2+eta.^4.*g.^4.*r.^4).*T.^4)+e.^4.*t.^4.*( ...
  1+(-1).*(1+3.*eta.^2.*g.^2.*r.^2).*T.^2+((-1)+2.*eta.^2.*g.^2.* ...
  r.^2+eta.^4.*g.^4.*r.^4).*T.^4+(1+eta.^2.*g.^2.*r.^2).*T.^6)+2.* ...
  e.^3.*t.^3.*(4+(-3).*(1+3.*eta.^2.*g.^2.*r.^2).*T.^2+2.*((-1)+2.* ...
  eta.^2.*g.^2.*r.^2+eta.^4.*g.^4.*r.^4).*T.^4+(1+eta.^2.*g.^2.* ...
 r.^2).*T.^6)+4.*eta.*g.*l.^7.*r.^2.*t.^(1/2).*(2+((-2)+e).*t).*T.* ...
  ((-4)+2.*((-2)+e).*t.*((-2)+((-2)+eta.^2.*g.^2.*r.^2).*T.^2)+((-2) ...
  +e).^2.*t.^2.*((-1)+((-2)+eta.^2.*g.^2.*r.^2).*T.^2+(3+eta.^2.* ...
  g.^2.*r.^2).*T.^4))+(-4).*eta.*g.*l.*r.^2.*t.^(1/2).*(2+e.*t).*T.* ...
  ((-4)+2.*e.*t.*((-2)+((-2)+eta.^2.*g.^2.*r.^2).*T.^2)+e.^2.*t.^2.* ...
  ((-1)+((-2)+eta.^2.*g.^2.*r.^2).*T.^2+(3+eta.^2.*g.^2.*r.^2).* ...
 T.^4))+l.^8.*(16+(-8).*((-2)+e).*t.*((-4)+(1+3.*eta.^2.*g.^2.* ...
r.^2).*T.^2)+4.*((-2)+e).^2.*t.^2.*(6+(-3).*(1+3.*eta.^2.*g.^2.* ...
 r.^2).*T.^2+((-1)+2.*eta.^2.*g.^2.*r.^2+eta.^4.*g.^4.*r.^4).*T.^4) ...
  +((-2)+e).^4.*t.^4.*(1+(-1).*(1+3.*eta.^2.*g.^2.*r.^2).*T.^2+((-1) ...
  +2.*eta.^2.*g.^2.*r.^2+eta.^4.*g.^4.*r.^4).*T.^4+(1+eta.^2.*g.^2.* ...
  r.^2).*T.^6)+2.*((-2)+e).^3.*t.^3.*(4+(-3).*(1+3.*eta.^2.*g.^2.* ...
  r.^2).*T.^2+2.*((-1)+2.*eta.^2.*g.^2.*r.^2+eta.^4.*g.^4.*r.^4).* ...
  T.^4+(1+eta.^2.*g.^2.*r.^2).*T.^6))+4.*eta.*g.*l.^3.*r.^2.*t.^( ...
  1/2).*T.*((-24)+3.*((-2)+e).*e.^2.*t.^3.*((-1)+((-2)+eta.^2.* ...
g.^2.*r.^2).*T.^2+(3+eta.^2.*g.^2.*r.^2).*T.^4)+4.*t.*(6+(-2).*(( ...
  -2)+r.^2+eta.^2.*g.^2.*r.^2).*T.^2+3.*e.*((-3)+((-2)+eta.^2.* ...
g.^2.*r.^2).*T.^2))+2.*e.*t.^2.*(3.*e.*((-3)+2.*((-2)+eta.^2.* ...
g.^2.*r.^2).*T.^2+(3+eta.^2.*g.^2.*r.^2).*T.^4)+(-2).*((-6)+((-8)+ ...
  r.^2+4.*eta.^2.*g.^2.*r.^2).*T.^2+(6+r.^2+2.*eta.^2.*g.^2.*r.^2).* ...
  T.^4)))+(-4).*eta.*g.*l.^5.*r.^2.*t.^(1/2).*T.*((-24)+3.*((-2)+e) ...
  .^2.*e.*t.^3.*((-1)+((-2)+eta.^2.*g.^2.*r.^2).*T.^2+(3+eta.^2.* ...
  g.^2.*r.^2).*T.^4)+4.*t.*(12+(-2).*((-4)+r.^2+2.*eta.^2.*g.^2.* ...
  r.^2).*T.^2+3.*e.*((-3)+((-2)+eta.^2.*g.^2.*r.^2).*T.^2))+2.* ...
t.^2.*(3.*e.^2.*((-3)+2.*((-2)+eta.^2.*g.^2.*r.^2).*T.^2+(3+ ...
eta.^2.*g.^2.*r.^2).*T.^4)+4.*((-3)+((-4)+r.^2+2.*eta.^2.*g.^2.* ...
 r.^2).*T.^2+(3+r.^2+eta.^2.*g.^2.*r.^2).*T.^4)+(-2).*e.*((-12)+(( ...
  -16)+r.^2+8.*eta.^2.*g.^2.*r.^2).*T.^2+(12+r.^2+4.*eta.^2.*g.^2.* ...
  r.^2).*T.^4)))+(-4).*l.^6.*(16+((-2)+e).^3.*e.*t.^4.*(1+(-1).*(1+ ...
  3.*eta.^2.*g.^2.*r.^2).*T.^2+((-1)+2.*eta.^2.*g.^2.*r.^2+eta.^4.* ...
  g.^4.*r.^4).*T.^4+(1+eta.^2.*g.^2.*r.^2).*T.^6)+(-4).*t.*(12+((-3) ...
  +(-3).*(1+3.*eta.^2.*g.^2).*r.^2+eta.^2.*g.^2.*r.^4).*T.^2+2.*e.*( ...
  (-4)+(1+3.*eta.^2.*g.^2.*r.^2).*T.^2))+4.*t.^2.*(e.^2.*(6+(-3).*( ...
  1+3.*eta.^2.*g.^2.*r.^2).*T.^2+((-1)+2.*eta.^2.*g.^2.*r.^2+ ...
eta.^4.*g.^4.*r.^4).*T.^4)+2.*(6+((-3)+(-3).*(1+3.*eta.^2.*g.^2).* ...
  r.^2+eta.^2.*g.^2.*r.^4).*T.^2+((-1)+(1+2.*eta.^2.*g.^2).*r.^2+ ...
  eta.^2.*g.^2.*(2+eta.^2.*g.^2).*r.^4).*T.^4)+(-1).*e.*(18+((-9)+( ...
  -3).*(1+9.*eta.^2.*g.^2).*r.^2+eta.^2.*g.^2.*r.^4).*T.^2+((-3)+(1+ ...
  6.*eta.^2.*g.^2).*r.^2+eta.^2.*g.^2.*(2+3.*eta.^2.*g.^2).*r.^4).* ...
  T.^4))+((-2)+e).^2.*t.^3.*((-4)+(3+(3+9.*eta.^2.*g.^2).*r.^2+(-1) ...
  .*eta.^2.*g.^2.*r.^4).*T.^2+(-2).*((-1)+(1+2.*eta.^2.*g.^2).*r.^2+ ...
  eta.^2.*g.^2.*(2+eta.^2.*g.^2).*r.^4).*T.^4+(-1).*(1+r.^2).*(1+ ...
  eta.^2.*g.^2.*r.^2).*T.^6+2.*e.*(4+(-3).*(1+3.*eta.^2.*g.^2.*r.^2) ...
  .*T.^2+2.*((-1)+2.*eta.^2.*g.^2.*r.^2+eta.^4.*g.^4.*r.^4).*T.^4+( ...
  1+eta.^2.*g.^2.*r.^2).*T.^6)))+(-4).*l.^2.*(16+((-2)+e).*e.^3.* ...
  t.^4.*(1+(-1).*(1+3.*eta.^2.*g.^2.*r.^2).*T.^2+((-1)+2.*eta.^2.* ...
  g.^2.*r.^2+eta.^4.*g.^4.*r.^4).*T.^4+(1+eta.^2.*g.^2.*r.^2).*T.^6) ...
  +(-4).*t.*(4+((-1)+(-3).*(1+eta.^2.*g.^2).*r.^2+eta.^2.*g.^2.* ...
 r.^4).*T.^2+2.*e.*((-4)+(1+3.*eta.^2.*g.^2.*r.^2).*T.^2))+4.*e.* ...
  t.^2.*((-6)+(3+(3+9.*eta.^2.*g.^2).*r.^2+(-1).*eta.^2.*g.^2.*r.^4) ...
  .*T.^2+(-1).*((-1)+(1+2.*eta.^2.*g.^2).*r.^2+eta.^2.*g.^2.*(2+ ...
 eta.^2.*g.^2).*r.^4).*T.^4+e.*(6+(-3).*(1+3.*eta.^2.*g.^2.*r.^2).* ...
  T.^2+((-1)+2.*eta.^2.*g.^2.*r.^2+eta.^4.*g.^4.*r.^4).*T.^4))+ ...
e.^2.*t.^3.*((-12)+(9+3.*(1+9.*eta.^2.*g.^2).*r.^2+(-1).*eta.^2.* ...
  g.^2.*r.^4).*T.^2+(-2).*((-3)+(1+6.*eta.^2.*g.^2).*r.^2+eta.^2.* ...
  g.^2.*(2+3.*eta.^2.*g.^2).*r.^4).*T.^4+(-1).*(3+r.^2).*(1+eta.^2.* ...
  g.^2.*r.^2).*T.^6+2.*e.*(4+(-3).*(1+3.*eta.^2.*g.^2.*r.^2).*T.^2+ ...
  2.*((-1)+2.*eta.^2.*g.^2.*r.^2+eta.^4.*g.^4.*r.^4).*T.^4+(1+ ...
eta.^2.*g.^2.*r.^2).*T.^6)))+2.*l.^4.*(48+3.*((-2)+e).^2.*e.^2.* ...
 t.^4.*(1+(-1).*(1+3.*eta.^2.*g.^2.*r.^2).*T.^2+((-1)+2.*eta.^2.* ...
  g.^2.*r.^2+eta.^4.*g.^4.*r.^4).*T.^4+(1+eta.^2.*g.^2.*r.^2).*T.^6) ...
  +(-8).*t.*(12+((-3)+(-3).*(2+3.*eta.^2.*g.^2).*r.^2+2.*eta.^2.* ...
  g.^2.*r.^4).*T.^2+3.*e.*((-4)+(1+3.*eta.^2.*g.^2.*r.^2).*T.^2))+ ...
  4.*t.^2.*(3.*e.^2.*(6+(-3).*(1+3.*eta.^2.*g.^2.*r.^2).*T.^2+((-1)+ ...
  2.*eta.^2.*g.^2.*r.^2+eta.^4.*g.^4.*r.^4).*T.^4)+(-2).*e.*(18+(( ...
  -9)+(-3).*(2+9.*eta.^2.*g.^2).*r.^2+2.*eta.^2.*g.^2.*r.^4).*T.^2+( ...
  (-3)+(2+6.*eta.^2.*g.^2).*r.^2+eta.^2.*g.^2.*(4+3.*eta.^2.*g.^2).* ...
  r.^4).*T.^4)+2.*(6+((-3)+(-3).*(2+3.*eta.^2.*g.^2).*r.^2+2.* ...
eta.^2.*g.^2.*r.^4).*T.^2+((-1)+2.*(1+eta.^2.*g.^2).*r.^2+(1+4.* ...
 eta.^2.*g.^2+eta.^4.*g.^4).*r.^4).*T.^4))+2.*e.*t.^3.*(3.*e.^2.*( ...
  4+(-3).*(1+3.*eta.^2.*g.^2.*r.^2).*T.^2+2.*((-1)+2.*eta.^2.*g.^2.* ...
  r.^2+eta.^4.*g.^4.*r.^4).*T.^4+(1+eta.^2.*g.^2.*r.^2).*T.^6)+2.*( ...
  12+((-9)+(-3).*(2+9.*eta.^2.*g.^2).*r.^2+2.*eta.^2.*g.^2.*r.^4).* ...
  T.^2+2.*((-3)+(2+6.*eta.^2.*g.^2).*r.^2+eta.^2.*g.^2.*(4+3.* ...
eta.^2.*g.^2).*r.^4).*T.^4+(3+2.*r.^2).*(1+eta.^2.*g.^2.*r.^2).* ...
 T.^6)+(-1).*e.*(36+((-27)+(-3).*(2+27.*eta.^2.*g.^2).*r.^2+2.* ...
eta.^2.*g.^2.*r.^4).*T.^2+2.*((-9)+2.*(1+9.*eta.^2.*g.^2).*r.^2+ ...
 eta.^2.*g.^2.*(4+9.*eta.^2.*g.^2).*r.^4).*T.^4+(9+2.*r.^2).*(1+ ...
 eta.^2.*g.^2.*r.^2).*T.^6)))))...
            ,0,inf,0,2*pi);
end

F = real(F);

end