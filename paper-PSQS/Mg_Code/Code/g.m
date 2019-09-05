function [ g ] = g( x )
    g = (x + 1) / 2 .* log2((x + 1) / 2) -...
        (x - 1) / 2 .* log2((x - 1) / 2);
g(isnan(g)) = 0;
end