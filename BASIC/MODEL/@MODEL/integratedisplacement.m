function [s,se]=integratedisplacement(S,q,varargin)

q=unfreevector(S,q);
[s,se]=integrate(S,@quadratureorder,@fun,q);

function I = fun(xi,elem,xnode,q)
N=calc_N(elem,xnode,xi);
qe=localize(elem,q);
I = N*qe;
return

function n=quadratureorder(elem)
n=orderN(elem);
return