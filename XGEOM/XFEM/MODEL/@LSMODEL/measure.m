function [s,se]=measure(S)

[s,se] = lsintegrate(S,0,@(xi,elem,xnode) 1,@(xi,elem,xnode) 0);


