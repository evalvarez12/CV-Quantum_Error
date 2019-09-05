function [Re_Val, Re_Dis, Ma_Dis] = Operation(T_Range, N, Para, A, z, TMin)
% ---- Iput ------
%
% T_Range: the range for the transmitivity T, when this parameter is set to
% 0, function will return the maximum distance.
% A: the average number of photon of the TMSV.
% B: channel noise.
% N: the number of the photon operated.
% th: the rate lower threshold.
% Tbs: the transmissivity for the beam splitter.
% Eln: flag for Eln calculation.
% sum_limit_n: summation limite for n.
%
% ---- Output ----
%
% Re_Val: return value, rate or Eln.
% Re_Dis: the corresponding distance of rates.
%
% ----------------
warning('off')
err = 0.1;
if nargin < 1
    Tbs = .9;
    T_Range = .5;
    A = 1.3;
    B = .001;
    N = 1;
    th = 1e-20;
    sum_limit_n = 30;
    sum_limit_m = 20;
    mode = 'CalRat';
    operation = 'RxS';
    TMin = 0;
    z = 0;
else
    %A = Para.A;
    B = Para.B;
    th = Para.th;
    Tbs = Para.Tbs;
    sum_limit_n = Para.sum_limit_n;
    sum_limit_m = Para.sum_limit_m;
    mode = Para.mode;
    operation = Para.operation;
end
%% Initialization
switch operation
    case 'TxA'
        Sumlet = zeros(sum_limit_n +1, sum_limit_n + N +1, sum_limit_m +1, sum_limit_m +1);
        r = zeros(sum_limit_n + N +1, sum_limit_n + N +1);
        rs = zeros(sum_limit_n + N +1, sum_limit_n + N +1);
    case 'RxA'
        Sumlet = zeros(sum_limit_n +1, sum_limit_n + N +1, sum_limit_m +1, sum_limit_m +1);
        r = zeros(sum_limit_n + N +1, sum_limit_n + N +1);
        rs = zeros(sum_limit_n + N +1, sum_limit_n + N +1);
    case 'TxS'
        sum_limit_n = sum_limit_n + N;
        Sumlet = zeros(sum_limit_n + 1, sum_limit_n + 1, sum_limit_m +1, sum_limit_m +1);
        r = zeros(sum_limit_n, sum_limit_n);
        rs = zeros(sum_limit_n, sum_limit_n);
    case 'RxS'
        sum_limit_n = sum_limit_n + N;
        Sumlet = zeros(sum_limit_n + 1, sum_limit_n + 1, sum_limit_m +1, sum_limit_m +1);
        r = zeros(sum_limit_n, sum_limit_n);
        rs = zeros(sum_limit_n, sum_limit_n);
    case 'No'
        Sumlet = zeros(sum_limit_n + 1, sum_limit_n + 1, sum_limit_m +1, sum_limit_m +1);
        r = zeros(sum_limit_n, sum_limit_n);
        rs = zeros(sum_limit_n, sum_limit_n);
end

% alpha, beta
n = 0:sum_limit_n+N;
if strcmp(mode, 'GetMaxA')
    A = 500;
    AL = 0.01;
    AR = 1000;
end

a = sqrt( ((A)/(1 + A)).^n / (1 + A));
b = sqrt( ((B)/(1 + B)).^n / (1 + B));

% z(n,k,m,l)
if z == 0
    z = zeros(sum_limit_n + 16, sum_limit_n + 16, sum_limit_m +1, sum_limit_m +1);
    for n = 0:sum_limit_n + 16
        for k = 0:n
            for m = 0:sum_limit_m
                for l = 0:m
                    z(n + 1, k + 1, m + 1, l + 1) = sqrt(nchoosek(n - k + l, l))...
                        * sqrt(nchoosek(k + m - l, k));
                end
            end
        end
    end
end

% rs(n,k)
for n = 0:sum_limit_n + sum_limit_m + N +1
    for k = 0:n
        rs(n + 1, k + 1) = sqrt(nchoosek(n, k))...
            * sqrt(Tbs)^(n - k) * sqrt(1 - Tbs)^k;
    end
end

% For transmitter case probability is irrelevant to transmissivity
P = 0;
switch operation
    case 'TxA'
        P = ((A+1)*(1-Tbs))^N/(1+A-A*Tbs)^(N+1);
    case 'TxS'
        P = ((A+1)*(1-Tbs))^N/(1+A-A*Tbs)^(N+1);
end

Ma_Dis = 0;
t = 1;
d = 1;
Distance = 0.0001;
T = 10^(Distance/-50);
delta = 1;
wentZero = 0;
calFlag = 0;
Dis_Res = 1e-3;
zPad = 1000;
if strcmp(mode, 'GetMinT')
    Re_Val = zeros(zPad,1);
end
%% Main Loop
while th
    switch mode
        case {'CalRat', 'CalEln'}
            th = 0;
            Re_Val = zeros(1, length(T_Range));
            Re_Dis = Re_Val;
    end
    for T_n = T_Range
        if ~strcmp(mode, 'GetMinT') % If not in the autonomous mode
            T = T_n;
        end
        %% Pre-alculation
%         B = T/(1-T) * err/2;
%         b = sqrt( ((B)/(1 + B)).^(0:sum_limit_n + N) / (1 + B));
        switch operation            
            case 'TxA'
                for n = 0:sum_limit_n + N
                    for k = 0:n
                        r(n + 1, k + 1) = sqrt(nchoosek(n, k))...
                            * sqrt(T)^(n - k) * sqrt(1 - T)^k;
                    end
                end
                
                for n = 0:sum_limit_n
                    for k = 0:n+N
                        for m = 0:sum_limit_m
                            for l = 0:m
                                Sumlet(n +1,k +1,m +1,l +1) = a(n +1) * rs(n+N +1, N +1)...
                                    * (-1)^k * r(n +N +1, k +1)...
                                    * b(m +1) * r(m +1, l +1) * z(n +N +1, k +1, m +1, l +1);
                            end
                        end
                    end
                end
        
                if strcmp(mode, 'GetMaxA')
                    % Re-calculate P
                    P = ((A+1)*(1-Tbs))^N/(1+A-A*Tbs)^(N+1);                    
                end
            case 'RxA'
                for n = 0:sum_limit_n
                    for k = 0:n
                        r(n + 1, k + 1) = sqrt(nchoosek(n, k))...
                            * sqrt(T)^(n - k) * sqrt(1 - T)^k;
                    end
                end
                for n = 0:sum_limit_n
                    for k = 0:n
                        for m = 0:sum_limit_m
                            for l = 0:m
                                Sumlet(n +1,k +1,m +1,l +1) = a(n +1) *rs(n-k+l+N +1, N +1)...
                                    * (-1)^k * r(n +1, k +1)...
                                    * b(m +1) * r(m +1, l +1) * z(n +1, k +1, m +1, l +1);
                            end
                        end
                    end
                end
            case 'TxS'
                for n = 0:sum_limit_n
                    for k = 0:n
                        r(n + 1, k + 1) = sqrt(nchoosek(n, k))...
                            * sqrt(T)^(n - k) * sqrt(1 - T)^k;
                    end
                end
                for n = N:sum_limit_n
                    for k = 0:n-N
                        for m = 0:sum_limit_m
                            for l = 0:m
                                Sumlet(n +1,k +1,m +1,l +1) = a(n +1)*(-1)^(k+N)*rs(n +1, N +1)*b(m +1)*r(n-N +1, k +1)...
                                    *r(m +1, l +1) * z(n -N +1,k +1,m +1,l +1);
                            end
                        end
                    end
                end
                % Re-calculate P
                if strcmp(mode, 'GetMaxA')
                    P = ((A)*(1-Tbs))^N/(1+A-A*Tbs)^(N+1);
                end
            case 'RxS'
                for n = 0:sum_limit_n
                    for k = 0:n
                        r(n + 1, k + 1) = sqrt(nchoosek(n, k))...
                            * sqrt(T)^(n - k) * sqrt(1 - T)^k;
                    end
                end
                for n = 0:sum_limit_n
                    for k = 0:n
                        for m = 0:sum_limit_m
                            for l = 0:m
                                if n-k+l-N>=0
                                    Sumlet(n +1,k +1,m +1,l +1) = a(n +1)*(-1)^(k+N)*r(n +1, k +1)*b(m +1)*r(m +1, l +1)...
                                        *z(n +1,k +1,m +1,l +1)*rs(n-k+l +1,N +1);
                                end
                            end
                        end
                    end
                end
            case 'No'
                for n = 0:sum_limit_n
                    for k = 0:n
                        r(n + 1, k + 1) = sqrt(nchoosek(n, k))...
                            * sqrt(T)^(n - k) * sqrt(1 - T)^k;
                    end
                end
                for n = 0:sum_limit_n
                    for k = 0:n
                        for m = 0:sum_limit_m
                            for l = 0:m
                                Sumlet(n +1,k +1,m +1,l +1) = ...
                                    a(n +1)*(-1)^k*r(n +1, k +1)*b(m +1)*r(m +1, l +1)...
                                    *z(n +1,k +1,m +1,l +1);
                            end
                        end
                    end
                end
        end
        %% Re_Val or Eln
        switch mode
            case 'CalEln' % Level of Entanglement Calculation Sequence
                switch operation
                    case 'TxA'
                        Left_Term = zeros(sum_limit_n + 1, sum_limit_n + sum_limit_m + N +1,...
                            sum_limit_n + N + sum_limit_m  +1, sum_limit_m +1);
                        for n = 0:sum_limit_n
                            for k = 0:n+N
                                for m = 0:sum_limit_m
                                    for l = 0:m
                                        Left_Term(n +1, n-k+l+N +1, k+m-l +1, m+1) =...
                                            Left_Term(n +1, n-k+l+N +1, k+m-l +1, m+1) +...
                                            Sumlet(n +1,k +1,m +1,l +1);
                                    end
                                end
                            end
                        end
                    case 'RxA'
                        Left_Term = zeros(sum_limit_n + 1, sum_limit_n + sum_limit_m + N + 1,...
                            sum_limit_n + sum_limit_m + 1, sum_limit_m + 1);
                        for n = 0:sum_limit_n
                            for k = 0:n
                                for m = 0:sum_limit_m
                                    for l = 0:m
                                        Left_Term(n +1, n-k+l+N +1, k+m-l +1, m+1) =...
                                            Left_Term(n +1, n-k+l+N +1, k+m-l +1, m+1) +...
                                            Sumlet(n +1,k +1,m +1,l +1);
                                    end
                                end
                            end
                        end
                        P.v = sum(Left_Term(Left_Term~=0).^2);
                    case 'TxS'
                        Left_Term = zeros(sum_limit_n + 1, sum_limit_n + sum_limit_m - N + 1,...
                            sum_limit_n + sum_limit_m + 1, sum_limit_m + 1);
                        for n = N:sum_limit_n
                            for k = 0:n-N
                                for m = 0:sum_limit_m
                                    for l = 0:m
                                        Left_Term(n +1, n-k+l-N +1, k+m-l +1, m+1) =...
                                            Left_Term(n +1, n-k+l-N +1, k+m-l +1, m+1) +...
                                            Sumlet(n +1,k +1,m +1,l +1);
                                    end
                                end
                            end
                        end
                    case 'RxS'
                        Left_Term = zeros(sum_limit_n + 1, sum_limit_n + sum_limit_m - N + 1,...
                            sum_limit_n + sum_limit_m + 1, sum_limit_m + 1);
                        for n = 0:sum_limit_n
                            for k = 0:n
                                for m = 0:sum_limit_m
                                    for l = 0:m
                                        if n-k+l-N>=0
                                            Left_Term(n +1, n-k+l-N +1, k+m-l +1, m+1) =...
                                                Left_Term(n +1, n-k+l-N +1, k+m-l +1, m+1) +...
                                                Sumlet(n +1,k +1,m +1,l +1);
                                        end
                                    end
                                end
                            end
                        end
                        P.v = sum(Left_Term(Left_Term~=0).^2);
                    case 'No'
                        Left_Term = zeros(sum_limit_n + 1, sum_limit_n + sum_limit_m + 1,...
                            sum_limit_n + sum_limit_m + 1, sum_limit_m + 1);
                        for n = 0:sum_limit_n
                            for k = 0:n
                                for m = 0:sum_limit_m
                                    for l = 0:m
                                        Left_Term(n +1, n-k+l +1, k+m-l +1, m+1) =...
                                            Left_Term(n +1, n-k+l +1, k+m-l +1, m+1)+...
                                            Sumlet(n +1,k +1,m +1,l +1);
                                    end
                                end
                            end
                        end
                        P.v = 1;
                end
                
                dim = size(Left_Term);
                sum_n = sum_limit_n;
                sum_k = dim(2);
                sum_l = sum_limit_m;
                sum_m = dim(3);
                
                q = zeros(sum_n * sum_k * sum_m * sum_l, 1);
                for n = 0:sum_n-1
                    for k = 0:sum_k-1
                        for m = 0:sum_m-1
                            for l = 0:sum_l-1
                                q(n*sum_k*sum_m*sum_l + k*sum_m*sum_l + m*sum_l + l + 1) = Left_Term(n +1,k +1,m +1,l +1)/sqrt(P.v);
                            end
                        end
                    end
                end
                
                pt = zeros(sum_n * sum_k, sum_n * sum_k);
                for i = 1:sum_m*sum_l
                    pt_left = q(i:sum_m*sum_l:end);
                    pt = pt + pt_left * pt_left';
                end
                
                ptt = zeros(sum_n * sum_k, sum_n * sum_k);
                for n = 1:sum_n
                    for m = 1:sum_n
                        ptt((-sum_k+1:0)+n * sum_k, (-sum_k+1:0)+m*sum_k)...
                            = pt((-sum_k+1:0)+n * sum_k, (-sum_k+1:0)+m*sum_k)';
                    end
                end
                eg = eig(ptt);
                Re_Val(t) = log2(1 - 2*sum(eg(eg<0)));
            case 'CalRat' % QKD Rate Calculation Sequence
                if T > TMin
                    Re_Val(t) = GetVar(Sumlet, sum_limit_n, sum_limit_m, N, P, operation, Tbs, A, B, T);
                end                
            case 'GetMinT' % Obtaining minimum transmissivity
                Re_Val(t) = GetVar(Sumlet, sum_limit_n, sum_limit_m, N, P, operation, Tbs, A, B, T);
                % Direction: 0 for increasing and 1 for decreasing
                if Re_Val(t) > 0
                    Re_Dis(d) = Distance;
                    d = d + 1;
                    % Here we want to speed up the approaching process
                    if ~wentZero
                        delta = 1 * delta;
                    end
                else
                    wentZero = 1;
                end
                
                % Approaching sequence
                if wentZero
                    if (Re_Val(t) <=0 || Re_Val(t) > th) % Stopping condition not satisfied
                        if Re_Val(t) * Re_Val(t-1) <= 0 % if changed symbol
                            delta = -delta/2;
                        end
                    else
                        Re_Val(end) = th;
                        th = 0;
                    end
                    
                    if abs(delta)<Dis_Res % OK enough resolution
                        Re_Dis(d) = Re_Dis(d-1);
                        Re_Val(end+1) = th;
                        th = 0;
                        Ma_Dis = Re_Dis(d);
                        break
                    end
                end
                
                Distance = Distance + delta;
                T = 10.^(-0.2 * Distance / 10);
            case 'GetMaxA' % Find the alpha that maximizes the rate given T sequence
                Rate = GetVar(Sumlet, sum_limit_n, sum_limit_m, N, P, operation, Tbs, A, B, T);
                Rate(Rate<=0) = 0;
                if ~calFlag
                    Ratemid = Rate;
                    Amid = A;
                    Amidmid = (Amid + AR)/2;
                    % Now have calculated Amid, go calculate Amidmid
                    A = Amidmid;
                    calFlag = 1;
                else
                    % Now have Amid and Amidmid, compare them
                    Ratemidmid = Rate;
                    if AL/AR > 0.995
                        Re_Val = Rate;
                        Re_Dis = A;
                        return
                    end
                    if Ratemidmid > Ratemid
                        AL = Amid;
                    else
                        AR = Amidmid;
                    end
                    Amid = (AL + AR)/2;
                    Amidmid = (Amid + AR)/2;
                    A = Amid;
                    calFlag = 0;
                end
                
                n = 0:sum_limit_n+N;
                a = sqrt( ((A)/(1 + A)).^n / (1 + A));
        end
        t = t + 1;
    end
end
Re_Val(Re_Val < 0) = 0;
% Re_Val = -Re_Val;
if strcmp(mode, 'GetMinT')
    Re_Val = Re_Val(Re_Val>0)';
    Re_Val = [Re_Val zeros(1, zPad - length(Re_Val))];
    Re_Dis = [Re_Dis zeros(1, zPad - length(Re_Dis))];
end
end

function [Rate, V_A] = GetVar(Sumlet, sum_limit_n, sum_limit_m, N, P, mode, Tbs, A, B, T)
V_A = 0;
V_B = 0;
V_E = 0;
V_F = 0;
Pt = 0;

C_AB = 0;
C_EF = 0;
C_EB = 0;
C_FB = 0;
switch mode
    case 'RxA' % RxA
        P = 0;
        for n = 0:sum_limit_n
            for k = 0:n
                for m = 0:sum_limit_m
                    for l = 0:m
                        for j = max(0, k-l):min(n, m+(k-l))
                            P = P + Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
                            V_A = V_A + n * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
                            V_B = V_B + (n+N-k+l) * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
                            V_E = V_E + (k+m-l) * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
                            V_F = V_F + m * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
                        end
                    end
                end
            end
        end
        
        for n = 0:sum_limit_n
            for k = 0:n
                for m = 1:sum_limit_m
                    for l = 0:m
                        for j = max(0, k-l):min(n, m+(k-l))
                            C_EF = C_EF + 2 * sqrt(m) * sqrt(k+m-l) * ...
                                Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m , j-(k-l) +1);
                        end
                    end
                end
            end
        end
        
        for n = 1:sum_limit_n
            for k = 0:n
                for m = 0:sum_limit_m
                    for l = 0:m
                        for j = max(0, k-l):min(n, m+(k-l))
                            C_AB = C_AB + 2 * sqrt(n) * sqrt(n+N-k+l) * ...
                                Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n, j +1, m +1, j-(k-l) +1);
                        end
                    end
                end
            end
        end
        

        C_FB = 0;
        for n = 0:sum_limit_n
            for k = 0:n
                for m = 1:sum_limit_m
                    for l = 0:m
                        for j = max(0, k-l+1):min(n, m+(k-l+1))
                            C_FB = C_FB + 2 * sqrt(n+N-k+l) * sqrt(m) * ...
                                Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m, j-(k-l+1) + 1);
                        end
                    end
                end
            end
        end
        
        C_EB = 0;
        for n = 0:sum_limit_n
            for k = 0:n
                for m = 0:sum_limit_m
                    for l = 0:m
                        for j = max(0, k-l+1):min(n, m+(k-l+1))
                            C_EB = C_EB + 2 * sqrt(n+N-k+l) * sqrt(m+k-l+1) * ...
                                Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l+1) + 1);
                        end
                    end
                end
            end
        end
%         C = A * T + (1 - T) * B;
%         Pt = (1 - Tbs)^N / (1 - C/(1+C) * Tbs)^(N + 1) / (1 + C);
%         P - Pt
%         V_Bt = 1 + 2 * (N * (C + 1) + C * Tbs) / (1 + C * (1 - Tbs));
%         V_At = 1 + 2 * (N + 1) * A/(1+A) * Tbs / (1 - C/(1+C) * Tbs) ...
%             * (1 - Tbs)^N / (1 - C/(1+C) * Tbs) ^ (N + 1) / (1 + C) / P
    case 'TxA' % TxA
%         for n = 0:sum_limit_n
%             for k = 0:n+N
%                 for m = 0:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l):min(n+N, m+(k-l))                            
%                             V_A = V_A + n * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                             V_B = V_B + (n+N-k+l) * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                             V_E = V_E + (k+m-l) * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                             V_F = V_F + m * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                         end
%                     end
%                 end
%             end
%         end
% 
%         for n = 0:sum_limit_n
%             for k = 0:n+N
%                 for m = 1:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l):min(n+N, m+(k-l))
%                             C_EF = C_EF + 2 * sqrt(m) * sqrt(k+m-l) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m , j-(k-l) +1);
%                         end
%                     end
%                 end
%             end
%         end
%         
%         for n = 1:sum_limit_n
%             for k = 0:n+N
%                 for m = 0:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l):min(n+N, m+(k-l))
%                             C_AB = C_AB + 2 * sqrt(n) * sqrt(n+N-k+l) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n, j +1, m +1, j-(k-l) +1);
%                         end
%                     end
%                 end
%             end
%         end
%         
%         C_FB = 0;
%         for n = 0:sum_limit_n
%             for k = 0:n+N
%                 for m = 1:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l+1):min(n+N, m+(k-l+1))
%                             C_FB = C_FB + 2 * sqrt(n+N-k+l) * sqrt(m) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m, j-(k-l+1) + 1);
%                         end
%                     end
%                 end
%             end
%         end
%         
%         C_EB = 0;
%         for n = 0:sum_limit_n
%             for k = 0:n+N
%                 for m = 0:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l+1):min(n+N, m+(k-l+1))
%                             C_EB = C_EB + 2 * sqrt(n+N-k+l) * sqrt(m+k-l+1) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l+1) + 1);
%                         end
%                     end
%                 end
%             end
%         end                
        
        V_A = 1 + 2 * (N + 1) * A * Tbs / (1 + A - A * Tbs);
        V_B1 = 1 + 2 * (N + A * (N + Tbs)) / (1 + A - A * Tbs);
        V_F = 1 + 2 * B;
        V_B = T * V_B1 + (1 - T) * V_F;
        V_E = (1 - T) * V_B1 + T * V_F;

        C_AB = sqrt(T) * 2 * sqrt(A*Tbs/(1+A)) * (N + 1) * (1 + A) / (1 + A - A * Tbs);
        C_EF = sqrt(T) * 2 * sqrt(B^2 + B);
        C_FB = sqrt(1 - T) * 2 * sqrt(B^2 + B);
        C_EB = sqrt(T*(1-T)) * (V_F - V_B1);
    case 'RxS' % RxS
        P = 0;
        for n = 0:sum_limit_n
            for k = 0:n
                for m = 0:sum_limit_m
                    for l = 0:m
                        if n-k+l-N>=0
                            for j = max(0, k-l):min(n, m+(k-l))
                                P = P + Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
                                V_A = V_A + n * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
                                V_B = V_B + (n-N-k+l) * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
                                V_E = V_E + (k+m-l) * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
                                V_F = V_F + m * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
                            end
                        end
                    end
                end
            end
        end
        
        for n = 0:sum_limit_n
            for k = 0:n
                for m = 1:sum_limit_m
                    for l = 0:m
                        if n-k+l-N>=0
                            for j = max(0, k-l):min(n, m+(k-l))
                                C_EF = C_EF + 2 * sqrt(m) * sqrt(k+m-l) * ...
                                    Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m , j-(k-l) +1);
                            end
                        end
                    end
                end
            end
        end
        
        for n = 1:sum_limit_n
            for k = 0:n
                for m = 0:sum_limit_m
                    for l = 0:m
                        if n-k+l-N>=0
                            for j = max(0, k-l):min(n, m+(k-l))
                                C_AB = C_AB + 2 * sqrt(n) * sqrt(n-N-k+l) * ...
                                    Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n, j +1, m +1, j-(k-l) +1);
                            end
                        end
                    end
                end
            end
        end
        
        C_FB = 0;
        for n = 0:sum_limit_n
            for k = 0:n
                for m = 1:sum_limit_m
                    for l = 0:m
                        if n-k+l-N>=0
                            for j = max(0, k-l+1):min(n, m+(k-l+1))
                                C_FB = C_FB + 2 * sqrt(n-N-k+l) * sqrt(m) * ...
                                    Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m, j-(k-l+1) + 1);
                            end
                        end
                    end
                end
            end
        end
        
        C_EB = 0;
        for n = 0:sum_limit_n
            for k = 0:n
                for m = 0:sum_limit_m
                    for l = 0:m
                        if n-k+l-N>=0
                            for j = max(0, k-l+1):min(n, m+(k-l+1))
                                C_EB = C_EB + 2 * sqrt(n-N-k+l) * sqrt(m+k-l+1) * ...
                                    Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l+1) + 1);
                            end
                        end
                    end
                end
            end
        end

%         C = T * A+ (1 - T) * B;
%         P = ((C)*(1-Tbs))^N/(1+C-C*Tbs)^(N+1); 
%         Tpc = C/(1+C)*Tbs;
%         Tpa = A/(1+A)*Tbs;
%         V_At = (N + Tpc)/(1-Tpc)^2*(C/(1+C)-Tpc)^N/(1-Tpc)^(N)/(1+C);
%         V_Att = (N*(A+1)+C*Tbs)*C^N*(1-Tbs)^N/(1+C-C*Tbs)^(N+2);
%         disp(['V_A = ' num2str(V_A)])
%         disp(['V_At = ' num2str(V_At)])
%         disp(['V_Att = ' num2str(V_Att)])       
        
        %V_Bt = (N + 1) * C * Tbs / (1 + C - C * Tbs)

    case 'TxS' % TxS
%         V_A = 0;
%         for n = N:sum_limit_n
%             for k = 0:n-N
%                 for m = 0:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l):min(n-N, m+(k-l))
%                             Pt = Pt + Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                             V_A = V_A + n * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                             V_B = V_B + (n-N-k+l) * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                             V_E = V_E + (k+m-l) * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                             V_F = V_F + m * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                         end
%                     end
%                 end
%             end
%         end
%  
%         for n = N:sum_limit_n
%             for k = 0:n-N
%                 for m = 1:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l):min(n, m+(k-l))
%                             C_EF = C_EF + 2 * sqrt(m) * sqrt(k+m-l) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m , j-(k-l) +1);
%                         end
%                     end
%                 end
%             end
%         end
%         
%         for n = N+1:sum_limit_n
%             for k = 0:n-N
%                 for m = 0:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l):min(n, m+(k-l))
%                             C_AB = C_AB + 2 * sqrt(n) * sqrt(n-N-k+l) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n, j +1, m +1, j-(k-l) +1);
%                         end
%                     end
%                 end
%             end
%         end
%         
%         C_FB = 0;
%         for n = N:sum_limit_n
%             for k = 0:n-N
%                 for m = 1:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l+1):min(n, m+(k-l+1))
%                             C_FB = C_FB + 2 * sqrt(n-N-k+l) * sqrt(m) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m, j-(k-l+1) + 1);
%                         end
%                     end
%                 end
%             end
%         end
% 
%         C_EB = 0;
%         for n = N:sum_limit_n
%             for k = 0:n-N
%                 for m = 0:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l+1):min(n, m+(k-l+1))
%                             C_EB = C_EB + 2 * sqrt(n-N-k+l) * sqrt(m+k-l+1) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l+1) + 1);
%                         end
%                     end
%                 end
%             end
%         end

        V_A = 1 + 2 * (N * A + N + A * Tbs)/(1 + A - A * Tbs);
        V_B1 = 1 + 2 * (N + 1) * A * Tbs / (1 + A - A * Tbs);
        V_F = 1 + 2 * B;
        V_B = T * V_B1 + (1 - T) * V_F;
        V_E = (1 - T) * V_B1 + T * V_F;

        C_AB = sqrt(T) * 2 * sqrt(A*Tbs/(1+A)) * (N + 1) * (1 + A)/(1 + A - A * Tbs);
        C_EF = sqrt(T) * 2 * sqrt(B^2 + B);
        C_FB = sqrt(1 - T) * 2 * sqrt(B^2 + B);
        C_EB = sqrt(T*(1-T)) * (V_F - V_B1);
    case 'No'  % No
        P = 1;
%         for n = 0:sum_limit_n
%             for k = 0:n
%                 for m = 0:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l):min(n, m+(k-l))
%                             V_A = V_A + Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                             V_B = V_B + (n-k+l) * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                             V_E = V_E + (k+m-l) * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                             V_F = V_F + m * Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l) +1);
%                         end
%                     end
%                 end
%             end
%         end
       
%                  
%         for n = 0:sum_limit_n
%             for k = 0:n
%                 for m = 1:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l):min(n, m+(k-l))
%                             C_EF = C_EF + 2 * sqrt(m) * sqrt(k+m-l) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m , j-(k-l) +1);
%                         end
%                     end
%                 end
%             end
%         end
%         
%         for n = 1:sum_limit_n
%             for k = 0:n
%                 for m = 0:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l):min(n, m+(k-l))
%                             C_AB = C_AB + 2 * sqrt(n) * sqrt(n-k+l) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n, j +1, m +1, j-(k-l) +1);
%                         end
%                     end
%                 end
%             end
%         end
%         
%         C_FB = 0;
%         for n = 0:sum_limit_n
%             for k = 0:n
%                 for m = 1:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l+1):min(n, m+(k-l+1))
%                             C_FB = C_FB + 2 * sqrt(n-k+l) * sqrt(m) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m, j-(k-l+1) + 1);
%                         end
%                     end
%                 end
%             end
%         end
%         
%         C_EB = 0;
%         for n = 0:sum_limit_n
%             for k = 0:n
%                 for m = 0:sum_limit_m
%                     for l = 0:m
%                         for j = max(0, k-l+1):min(n, m+(k-l+1))
%                             C_EB = C_EB + 2 * sqrt(n-k+l) * sqrt(m+k-l+1) * ...
%                                 Sumlet(n +1, k +1, m +1, l +1) * Sumlet(n +1, j +1, m +1, j-(k-l+1) + 1);
%                         end
%                     end
%                 end
%             end
%         end
        %disp(['1: ' num2str(C_EB)]);
        V_A = 2 * A + 1;
        V_F = 2 * B + 1;
        V_B = T * V_A + (1 - T) * V_F;
        V_E = (1 - T) * V_A + T * V_F;
        
        C_AB = sqrt(T) * 2 * sqrt(A^2 + A);
        C_EF = sqrt(T) * 2 * sqrt(B^2 + B);
        C_FB = sqrt(1 - T) * 2 * sqrt(B^2 + B);
        C_EB = sqrt(T*(1-T)) * (V_F - V_A);
        %disp(['2: ' num2str(C_EB)])
end

switch mode
    case {'RxA', 'RxS'}
        
        V_A = 1 + 2 * V_A / P;
        V_B = 1 + 2 * V_B / P;
        V_F = 1 + 2 * V_F / P;
        V_E = 1 + 2 * V_E / P;
        
        C_AB = C_AB / P;
        C_EF = C_EF / P;
        C_FB = C_FB / P;
        C_EB = C_EB / P;

end

%% Entropy and rate calculation
[Rate] =  CM2Rate(V_A, V_B, C_AB, V_E, V_F, C_EF, C_FB, C_EB);
Rate = P * Rate;
end
