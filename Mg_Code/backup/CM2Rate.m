function [Rate, I_AB, X_AgB] = CM2Rate(V_A, V_B, C_AB)
% Reconciliation efficiency
f = 0.95;
% Detection efficiency
eta = 1;
% Var(thermal)
nu = 0;
% Detection effifiency and thermal noise
V_Bn = eta * V_B + (1-eta) * (1 + 2 * nu);
C_ABn = sqrt(eta) * C_AB;

% S(AB) - S(A|B), trusted device scenario
E_AB = .5 * sqrt((V_A + V_B)^2 - 4 * C_AB^2) + .5 * [V_B-V_A V_A-V_B];
E_AgB = sqrt(V_A*(V_A - C_ABn^2 / (V_Bn)));
X_AgB = sum(g(E_AB)) - sum(g(E_AgB));

V_BgA = V_Bn - C_ABn ^ 2 / (V_A);
I_AB = f * 0.5 * log2(V_Bn/V_BgA);
Rate = I_AB - abs(X_AgB);
end

% I = diag([1,1]);
% S = diag([1,-1]);
% Pi = diag([1,0]);
% O_u = [0,1;-1,0];
% O_n = [0,0;0,0];
% O = [O_u,O_n;O_n,O_u];

% S(EF) - S(EF|B)
% E_EF = .5 * sqrt((V_E + V_F)^2 - 4 * C_EF^2) + .5 * [V_E-V_F V_F-V_E];
% L_EFB = [C_EB * I; C_FB * S];
% CM_EFgB = [V_E * I, C_EF * S;C_EF * S, V_F * I] - L_EFB * Pi/V_B * L_EFB';
% E_EFgB = abs(eig(1i*O*CM_EFgB));
% X_EFgB = P * real(sum(g(E_EF)) - sum(g(E_EFgB))/2);
% disp(['S(EF) - S(EF|B) = ' num2str(X_EFgB)])
 
% % S(EF) - S(E|FB)
% E_EgF = sqrt(V_E * (V_E - C_EF^2 / (V_F)));
% V_EgB = V_E - C_EB^2/V_B;
% V_FgB = V_F - C_FB^2/V_B;
% C_EFgB = C_EF - C_EB*C_FB/V_B;
% V_EgFB = V_EgB - C_EFgB^2 /V_FgB;
% E_EgFB = sqrt((V_E - C_EB^2/V_B) * (V_E - C_EF^2/V_F));
% X_EgFB = P * (sum(g(E_EF)) - g(E_EgFB));
% disp(['S(EF) - S(E|FB) = ' num2str(X_EgFB)])

% Rate
% V_BgA = V_B - C_AB ^ 2 / (V_A);
% I_AB = f * 0.5 * log2(V_B/V_BgA);