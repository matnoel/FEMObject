function R=funV(R,fun,varargin)
fun = fcnchk(fun);

for k=1:R.m
   R.V{k}=fun(R.V{k},varargin{:});
end
