function f = getparam(u,field)
% function f = getparam(u,field)

if nargin==1
    f = u.param;
else
    f = getfield(u.param,field); 
end
