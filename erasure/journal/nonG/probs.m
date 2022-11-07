% % Coefficients (with 95% confidence bounds):
% a1 =     0.07544;
% b1 =       4.059;
% c1 =       1.178;
% a2 =     0.04125;
% b2 =       4.466;
% c2 =       2.569;
% a3 =     -0.0105;
% b3 =       5.892;
% c3 =       1.982;
% a4 =     0.01371;
% b4 =       3.915;
% c4 =       5.479;
% a5 =     0.01022;
% b5 =       4.643;
% c5 =       18.87;
% 
% b1 = b1 +20;
% b2 = b2 +20;
% b3 = b3 +20;
% b4 = b4 +20;
% b5 = b5 +20;
% 
% % Total area of the function
% N = sqrt(pi)*(a1*c1 + a2*c2 + a3*c3 + a4*c4 + a5*c5);
% 
% % Negative side
% Pn = sqrt(pi)/2 * ((a1*c1*(-erf(b1/c1) + 1)) + (a2*c2*(-erf(b2/c2) + 1)) + (a3*c3*(-erf(b3/c3) + 1)) + ...
%     (a4*c4*(-erf(b4/c4) + 1)) +(a5*c5*(-erf(b5/c5) + 1)));
% 
% % Positive side
% Pp = sqrt(pi)/2 * ((a1*c1*(erf(b1/c1) + 1)) + (a2*c2*(erf(b2/c2) + 1)) + (a3*c3*(erf(b3/c3) + 1)) + ...
%     (a4*c4*(erf(b4/c4) + 1)) +(a5*c5*(erf(b5/c5) + 1)));
% 
% 
% Pn^2
% (Pn/N)^2

x = gg(:,: ,1);
y = gg(:,:,2);

surf(x,y,ps);
% contour(x,y,p);
N = sum(ps, 'all');
Ny = sum(ps);
sum(Ny(17:end))/N