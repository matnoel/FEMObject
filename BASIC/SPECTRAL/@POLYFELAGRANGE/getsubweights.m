function x = getsubweights(h,i)
% function x = getsubweights(h,i)
x=getparam(h,'subweights');
if nargin==2
    x = x{i};
end


