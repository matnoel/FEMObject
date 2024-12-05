function plot(ls,M,varargin)
% function plot(ls,M)

if nargin==1 || ~isa(M,'MODEL')
    error('rentrer un MODEL en deuxieme argument')
end

if ~iseval(ls)
    ls = lseval(ls,M);
end

if numel(ls)>1
    error('ls est une multilevelset')
end

plot(FENODEFIELD(ls.value),M,varargin{:})
