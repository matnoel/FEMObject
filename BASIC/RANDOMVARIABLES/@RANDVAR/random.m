function A = random(u,varargin)
% function A = random(u,varargin)

n = getclassin('double',varargin,{1});
if ~isa(n,'cell')
    n = {n};
end

if ischarin('init',varargin);
    initstate
end

A = random(u.type,u.param{:,2},n{:});
