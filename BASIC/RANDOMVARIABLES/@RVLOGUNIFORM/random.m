function A = random(u,varargin)
% function A = random(u,varargin)

[rep,pos] = isclassin('double',varargin);
rep = ischarin('init',varargin);
if rep
    initstate
end
param = getparam(u);

a = log(param.A);
b = log(param.B);

A = exp(random('uniform',a,b,varargin{pos}));
