function [a,xi] = random(x,n1,n2)
% function [a,xi] = random(x,n1,n2)

if nargin<=1
    n1 = 1;
end
if nargin<=2
    n2 = 1;
end
n = [n1,n2];
xi = random(RANDVARS(x),prod(n),1);
xi = [xi{:}];
a = randomeval(x,xi);
