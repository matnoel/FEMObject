function x = getreppoints(h,i)
% function x = getreppoints(h,i)
x=getparam(h,'reppoints');
if nargin==2
    x = x{i};
end


