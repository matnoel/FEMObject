function param = funallparam(param,fun,varargin)
% function param = funallparam(param,fun,varargin)

fun = fcnchk(fun);
for i=1:size(param,1)
    param{i,2} = fun(param{i,2},varargin{:});
end


