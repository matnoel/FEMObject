function [a,xi] = random(apc,varargin)
% function [a,xi] = random(apc,varargin)

rep = ischarin('init',varargin);
if rep
    initstate
end
n = getclassin('double',varargin);
if isa(n,'cell')
    n = [n{:}];
elseif isempty(n)
    n = 1;
end

[D,xi] = random(apc.POLYCHAOS,prod(n));

if iscell(apc)
    a = getvalue(multimtimes(D,apc.MULTIMATRIX));
else
    a = double(apc) * D';
end

if n==1
    if iscell(apc)
        a = a{1};
    end
    a = reshape(a,size(apc));
elseif all(size(apc)==1)
    if length(n)==1
        n = [1,n];
    end
    if iscell(apc)
        a = [a{:}];
    end
    a = reshape(a,n);
else
    if length(n)==1
        n = [n,1];
    end
    a = MULTIMATRIX(a,size(apc),n);
end
