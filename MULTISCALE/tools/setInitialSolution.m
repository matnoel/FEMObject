function u0 = setInitialSolution(initType,b,u_old,varargin)
% function u = setInitialSolution(initType,b,u_old,varargin)
% Sets initial solution according to initType
% initType: type of initialization ('zero', 'one' or 'previter'), 'zero' by default

if isnumeric(initType)
    u0 = initType;
else
    switch lower(initType)
        case 'zero'
            u0 = zeros(size(b));
        case 'one'
            u0 = ones(size(b));
        case 'previter'
            u0 = u_old;
        otherwise
            error('Wrong argument initType')
    end
end

end
