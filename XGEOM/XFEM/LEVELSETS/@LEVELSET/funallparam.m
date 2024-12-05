function ls = funallparam(ls,fun,varargin)
% function ls = funparam(ls,fun,varargin)
% appliquer la fonction fun aux parametres de la levelset
% varargin sont les arguments de fun

fun = fcnchk(fun);
for i=1:size(ls.value,1)
    ls.value{i,2} = fun(ls.value{i,2},varargin{:});
end

