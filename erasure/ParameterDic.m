classdef ParameterDic
    properties
        Size
        Pars1
        Pars2
        Pars3
        Pars12
        Pars13
        Pars23
        Vs
    end

    methods
        function obj=setupImpl(obj)
            file_F = load('data/dic_Fs');
            file_pars = load('data/dic_Fs_pars2');
            all_pars = file_pars.results_pars2;
            fs = file_F.results;
            [l, ~] = size(fs);
            obj.Size = l;
            obj.Vs = fs(:, 1);
            obj.Pars1 = all_pars(:, 1);
            obj.Pars2 = all_pars(:, 2);
            obj.Pars3 = all_pars(:, 3);
            obj.Pars12 = all_pars(:, 4);
            obj.Pars13 = all_pars(:, 5);
            obj.Pars23 = all_pars(:, 6);
        end
        
        
        
        function [F] = Fidelity(obj, V, modes, sigma)
            [~, ind] = min(abs(obj.Vs - V));
%             disp([ind]);
            switch modes
                case '1'
                    F = fid_tmsv_ri(V, obj.Pars1(ind), '1', sigma);
%                     disp([obj.Pars1(ind, 1) obj.Pars1(ind, 2)]);
                case '2'
                    F = fid_tmsv_ri(V, obj.Pars2(ind), '2', sigma);
%                     disp([obj.Pars2(ind, 1) obj.Pars2(ind, 2)]);
                case '3'
                    F = fid_tmsv_ri(V, obj.Pars3(ind), '3', sigma);
%                     disp([obj.Pars3(ind, 1) obj.Pars3(ind, 2)]);
                case '12'
                    F = fid_tmsv_ri(V, obj.Pars12(ind), '12', sigma);
%                     disp([obj.Pars12(ind, 1) obj.Pars12(ind, 2)]);
                case '13'
                    F = fid_tmsv_ri(V, obj.Pars13(ind), '13', sigma);
%                     disp([obj.Pars13(ind, 1) obj.Pars13(ind, 2)]);
                case '23'
                    F = fid_tmsv_ri(V, obj.Pars23(ind), '23', sigma);
%                     disp([obj.Pars23(ind, 1) obj.Pars23(ind, 2)]);
            end
        end
        
        
        function [F] = Full_fidelity(obj, Pe, V, sigma)
            F_1 = obj.Fidelity(V, '1', sigma);
            F_2 = obj.Fidelity(V, '2', sigma);
            F_3 = obj.Fidelity(V, '3', sigma);
            F_12 = obj.Fidelity(V, '12', sigma);
            F_13 = obj.Fidelity(V, '13', sigma);
            F_23 = obj.Fidelity(V, '23', sigma);
            F_123 = fid_tmsv(1, 1, 1, '123', sigma);
            
            F = (1-Pe)^3  + Pe*(1-Pe)^2*F_1 + + Pe*(1-Pe)^2*F_2 + + Pe*(1-Pe)^2*F_3 + ...
                Pe^2*(1-Pe)*F_12 + Pe^2*(1-Pe)*F_13 +Pe^2*(1-Pe)*F_23 + Pe^3*F_123;
        end
        
        function [F] = Partial_fidelity(obj, Pe, V, sigma)
            F_1 = obj.Fidelity(V, '1', sigma);
            F_2 = obj.Fidelity(V, '2', sigma);
            F_3 = obj.Fidelity(V, '3', sigma);

            
            F = (1-Pe)^3  + Pe*(1-Pe)^2*F_1 + + Pe*(1-Pe)^2*F_2 + + Pe*(1-Pe)^2*F_3;
        end
    end
end
        
        