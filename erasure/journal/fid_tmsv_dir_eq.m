function [F] = fid_tmsv_dir_eq(T, eps, sigma)

F =  2 ./ (2*sigma*(1-sqrt(T)).^2 + 2 + eps);
F = real(F);
end