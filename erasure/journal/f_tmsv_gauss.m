function [F] = f_tmsv_gauss(a1, b1, c1, a2, b2, c2)
    V1 = [a1 0 c1 0 ; 0 a1 0 -c1 ; c1 0 b1 0 ; 0 -c1 0 b1];
    V2 = [a2 0 c2 0 ; 0 a2 0 -c2 ; c2 0 b2 0 ; 0 -c2 0 b2];
    O = [0 1 0 0 ; -1 0 0 0 ; 0 0 0 1 ; 0 0 -1 0];

    Vaux = transpose(O) * inv(V1 + V2) * (O/4 + V2*O*V1);

    Vi = inv(Vaux*O);

    F4tot = det(2*(sqrtm(eye(4) + (Vi*Vi)/4) + eye(4))*Vaux);

    F = nthroot(F4tot/det(V1+V2),4);

end