function A = random(u,varargin)
% function A = random(u,varargin)

[rep,pos] = isclassin('double',varargin);
rep = ischarin('init',varargin);
if rep
    initstate
end
param = getparam(u);

A = param.x0+(param.x1-param.x0)*random('beta',param.a,param.b,varargin{pos});
