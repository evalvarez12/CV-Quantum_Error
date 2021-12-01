clc;clear



num = 12;
Ps = linspace(1, .6, num);

F = zeros(1, num);
Fdir = F;


sigma = 10;
F123 = fid_tmsv(1, 1, 1, '123', sigma);


parfor i = 1:num
    disp([i]);
    
    F(i) = full_opt_fid_tmsv(Pe);
    
    Fdir(i) = (1-Pe) + Pe * F123;
    

end

results = [F(:), Fdir(:)];
save('results_fullopt_Ftmsv.mat', 'results');

