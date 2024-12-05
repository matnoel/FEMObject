function param=funrandomparam(param,fun,varargin)

fun = fcnchk(fun);
for i=1:size(param,1)
    if israndom(param{i,2})
        param{i,2} = fun(param{i,2},varargin{:});
    end
end
