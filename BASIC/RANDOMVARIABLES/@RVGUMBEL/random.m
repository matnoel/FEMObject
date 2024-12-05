function A = random(u,varargin)
% function A = random(u,varargin)

[rep,pos] = isclassin('double',varargin);
rep = ischarin('init',varargin);
if rep
    initstate
end
param = struct2cell(getparam(u));

A = gumbelrnd(param{:},varargin{pos});
