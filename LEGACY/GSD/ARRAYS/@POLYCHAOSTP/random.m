function [D,xi] = random(PC,varargin)
% function [D,xi] = random(PC,varargin)

[rep,pos] = isclassin('double',varargin);
if rep
    n = [varargin{pos}];
else
    n = 1;    
end
n = prod(n);
rv = RANDVARS(PC);
xi = random(rv,n,1);
xi = [xi{:}];
D = polyval(PC,xi);
