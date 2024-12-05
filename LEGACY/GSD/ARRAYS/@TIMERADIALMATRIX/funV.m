function R=funV(R,fun,varargin)
fun = fcnchk(fun);

R.V = fun(R.V,varargin{:});
