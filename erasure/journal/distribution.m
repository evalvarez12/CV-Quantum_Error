function [ts] = distribution(pd, T0, mu, sigma, N)
    s = rand(N,1);
    s(s>pd) = 1;
    s(s ~= 1) = T0;
    
    
    ts = normrnd(mu, sigma, N, 1);
    ts = ts .* s;
    
    ts = max(ts, 0);
    
    
    extra = length(ts(ts>1));
    i = randi([1 N-extra], extra, 1);
    support = ts(ts <= 1);
    
    ts(ts > 1) = support(i);
        
%     histogram(ts, 100);
end