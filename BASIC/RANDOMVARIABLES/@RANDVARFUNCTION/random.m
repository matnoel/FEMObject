function x = random(rv,varargin)
% function x = random(rv,varargin)

p = random(rv.RV,varargin{:});
rv.param(rv.randomparam) = p;

x = rv.fun(rv.param{:});
