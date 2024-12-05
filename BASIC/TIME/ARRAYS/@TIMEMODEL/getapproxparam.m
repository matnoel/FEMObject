function p = getapproxparam(T,r)
% function p = getapproxparam(T,r)

if nargin==2
    p = getparam(T.approxparam,r);
else
    p = T.approxparam;
end

