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
            file_pars = load('data/dic_Fs_pars');
            all_pars = file_pars.results_pars;
            fs = file_F.results;
            [l, ~] = size(fs);
            obj.Size = l;
            obj.Vs = fs(:, 1);
            obj.Pars1 = [all_pars(:, 1) all_pars(:, 2)];
            obj.Pars2 = [all_pars(:, 3) all_pars(:, 4)];
            obj.Pars3 = [all_pars(:, 5) all_pars(:, 6)];
            obj.Pars12 = [all_pars(:, 7) all_pars(:, 8)];
            obj.Pars13 = [all_pars(:, 9) all_pars(:, 10)];
            obj.Pars23 = [all_pars(:, 11) all_pars(:, 12)];
        end
        
        
        
        function [F] = Fidelity(obj, V, modes, sigma)
            [~, ind] = min(abs(obj.Vs - V));
%             disp([ind]);
            switch modes
                case '1'
                    F = fid_tmsv(V, obj.Pars1(ind, 1), obj.Pars1(ind, 2), '1', sigma);
%                     disp([obj.Pars1(ind, 1) obj.Pars1(ind, 2)]);
                case '2'
                    F = fid_tmsv(V, obj.Pars2(ind, 1), obj.Pars2(ind, 2), '2', sigma);
%                     disp([obj.Pars2(ind, 1) obj.Pars2(ind, 2)]);
                case '3'
                    F = fid_tmsv(V, obj.Pars3(ind, 1), obj.Pars3(ind, 2), '3', sigma);
%                     disp([obj.Pars3(ind, 1) obj.Pars3(ind, 2)]);
                case '12'
                    F = fid_tmsv(V, obj.Pars12(ind, 1), obj.Pars12(ind, 2), '12', sigma);
%                     disp([obj.Pars12(ind, 1) obj.Pars12(ind, 2)]);
                case '13'
                    F = fid_tmsv(V, obj.Pars13(ind, 1), obj.Pars13(ind, 2), '13', sigma);
%                     disp([obj.Pars13(ind, 1) obj.Pars13(ind, 2)]);
                case '23'
                    F = fid_tmsv(V, obj.Pars23(ind, 1), obj.Pars23(ind, 2), '23', sigma);
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
    end
end
        
        