function [s,se]=integrateepsilon(S,q,varargin)

q=unfreevector(S,q);
[s,se]=integrate(S,@quadratureorder,@fun,q);

function I = fun(xi,elem,xnode,q)
B=calc_B(elem,xnode,xi);
qe=localize(elem,q);
I = B*qe;
return
 
function n=quadratureorder(elem)
n=orderB(elem);
return