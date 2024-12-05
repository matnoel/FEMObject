function [se,fieldtype,fieldstorage,fieldddl] = calc_elemfields(S,fun,varargin)
% function [se,fieldtype,fieldstorage,fieldddl] = calc_elemfields(S,fun,varargin)
% fun : methode de la classe ELEMENT , varargin : ses arguments
% function [se,fieldtype,fieldstorage,fieldddl] = calc_elemfields(S,funname,varargin)
% funname = 'sigma', 'epsilon', 'energyint', 'fint'

display_ = ischarin('display',varargin);
if display_
    fprintf('\n COMPUTING FINITE ELEMENT FIELD \n')
end

if isa(fun,'char')
    fun = eval(['@' fun]);
end
fun = fcnchk(fun);

for p=1:S.nbgroupelem
    if display_
        fprintf('-> Computing element field of element group %3d / %3d ... ',p,S.nbgroupelem)
    end
    [temp,fieldtype,fieldstorage,fieldddl] = fun(S.groupelem{p},getnode(S),varargin{:});
    for i=1:length(temp)
        se{i}{p} = temp{i};
    end
    if display_
        fprintf('\n')
    end
end

for i=1:length(se)
    se{i} = FEELEMFIELD(se{i},'storage',fieldstorage,'type',fieldtype,'ddl',fieldddl);
end