function [se,fieldtype,fieldstorage]=calc_elemfield(S,fun,varargin)
% function M=calc_elemfield(S,fun,varargin)
% fun :  methode de la classe ELEMENT , varargin : ses arguemnts
% function M=calc_elemfield(S,funname,varargin)
% funname = 'sigma' 'epsilon'

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
    [se{p},fieldtype,fieldstorage,fieldddl] = fun(S.groupelem{p},getnode(S),varargin{:});
    if display_
        fprintf('\n')
    end
end

se = FEELEMFIELD(se,'storage',fieldstorage,'type',fieldtype,'ddl',fieldddl);
