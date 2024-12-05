function [s,se]=integrateenergyproduct(S,q1,q2,varargin)

q1=unfreevector(S,q1);
q2=unfreevector(S,q2);
[s,se]=integrate(S,@quadratureorder,@fun,q1,q2,varargin{:});

function I = fun(xi,elem,xnode,q1,q2,varargin)

mat = getmaterial(elem);
C = calc_opmat(mat,elem,xnode,xi);
B=calc_B(elem,xnode,xi);
qe1=localize(elem,q1);
qe2=localize(elem,q2);
ep1 = B*qe1;
ep2 = B*qe2;
I = ep1'*C*ep2;
return
 
function n=quadratureorder(elem)
n=2*orderB(elem);
return