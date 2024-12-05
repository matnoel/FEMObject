function x = getsubpoints(h,i)
% function x = getsubpoints(h,i)
x=getparam(h,'subpoints');
if nargin==2
    x = x{i};
end


