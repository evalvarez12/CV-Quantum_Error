function [x]  = gfun(v)
    if v == 1
        x = 0;
    else
        x=(v+1)/2.*log2((v+1)/2) - (v-1)/2.*log2((v-1)/2);
    end
end