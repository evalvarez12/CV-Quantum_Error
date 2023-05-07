function [res]=RCI_dir(vs, T, eps, pe)
    
    a = vs;
    b = T.*(vs-1) + eps +1;
    c = sqrt(T).*sqrt(vs^2-1);

    z = sqrt((a+b).^2 - 4*c.^2);
    v1 = (z+(a-b))/2;
    v2 = (z-(a-b))/2;
    
    r  = acosh(vs)/2;
    nth = sinh(r)^2;


    Sth = gfun(vs);
%     Sth = gfun(2*nth +1);
    S = real(gfun(v1) + gfun(v2));
    
    S(S<0) = 0;
    
%     v1 
%     v2
% 
%     Sth
%     S

    res=(1-pe)*Sth  - (1-pe)*S;

end