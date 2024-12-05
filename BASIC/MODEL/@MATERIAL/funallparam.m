function mat = funallparam(mat,fun,varargin)
% function ls = funparam(ls,fun,varargin)
% appliquer la fonction fun aux parametres du materiau
% varargin sont les arguments de fun

fun = fcnchk(fun);
for i=1:size(mat.param,1)
    mat.param{i,2} = fun(mat.param{i,2},varargin{:});
end

