function x = getweights(h,i)
% function x = getweights(h,i)
x=getparam(h,'weights');
if nargin==2
    x = x{i};
end


