function p = getevolparam(T,r)
% function p = getevolparam(T,r)

if nargin==2
    p = getparam(T.evolparam,r);
else
    p = T.evolparam;
end

