function [D,xi] = lhsrandom(PC,varargin)
% function [D,xi] = lhsrandom(PC,varargin)

[rep,pos] = isclassin('double',varargin);
if rep
    n = [varargin{pos}];
else
    n = 1;    
end
n = prod(n);
rv = RANDVARS(PC);
xi = lhsrandom(rv,n,1);
xi = [xi{:}];
D = polyval(PC,xi);
