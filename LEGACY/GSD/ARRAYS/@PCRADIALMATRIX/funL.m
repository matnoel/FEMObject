function R=funL(R,fun,varargin)
fun = fcnchk(fun);

R.L = fun(R.L,varargin{:});
%for k=1:R.m
%   R.L{k}=fun(R.L(k),varargin{:});
%end