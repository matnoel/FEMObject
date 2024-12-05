function [s,se]=measure(S)

[s,se] = integrate(S,0,@(xi,elem,xnode) 1);


