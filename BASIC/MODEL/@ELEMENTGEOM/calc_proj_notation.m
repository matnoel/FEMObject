function P = calc_proj_notation(elem)
% function P = calc_proj_notation(elem)
% Change of notation matrix or transition matrix P 
% for strain vector epsilon from Kelvin-Mandel notation to Voigt notation 
% and for stress vector sigma from Voigt notation to Kelvin-Mandel notation 
% epsilon_V = P * epsilon_K and sigma_K = P * sigma_V 

dim = getdim(elem);

switch dim
    case 1
        P = 1;
    case 2
        if isaxi(elem)
            P = diag([1,1,1,sqrt(2)]);
        else
            P = diag([1,1,sqrt(2)]);
        end
    case 3
        P = diag([1,1,1,sqrt(2),sqrt(2),sqrt(2)]);
end

end
