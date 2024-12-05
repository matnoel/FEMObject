function [s,se]=integrateelasticenergy(S,q,varargin)

q=unfreevector(S,q);
[s,se]=integrate(S,@quadratureorder,@fun,q);

function I = fun(xi,elem,xnode,q)
mat = getmaterial(elem);
C = calc_opmat(mat,elem);
B=calc_B(elem,xnode,xi);
qe=localize(elem,q);
ep = B*qe;
I = 1/2*ep'*C*ep;
return
 
function n=quadratureorder(elem)
n=2*orderB(elem);
return