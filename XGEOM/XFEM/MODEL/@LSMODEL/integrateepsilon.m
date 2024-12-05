function [s,se]=integrateepsilon(S,q,varargin)

q=unfreevector(S,q);
[s,se]=lsintegrate(S,@quadratureorder,@fun,@funout,q);

function I = fun(xi,elem,xnode,q)
B=calc_Bls(elem,xnode,xi);
qe=localize(elem,q);
I = B*qe;
return

function I = funout(xi,elem,xnode,q)
I = zeros(getnbddlpergauss(elem),1);
return
 
 
function n=quadratureorder(elem)
n=orderB(elem);
return