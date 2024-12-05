function param=funparam(param,num,fun,varargin)
% function param=funparam(param,num,fun,varargin)

fun = fcnchk(fun);
if isa(num,'char')  || isa(num,'cell')
    [rep,num] = ischarin(num,param(:,1));
    num = num(find(rep));
end

for i=1:length(num)
    param{num(i),2} = fun(param{num(i),2},varargin{:});
end


