function [s,se]=inertia_matrix(S,xG)
if nargin==1
    xG=zeros(1,getindim(S));
end
[s,se]=integrate(S,@quadratureorder,@fun,xG);

function I = fun(xi,elem,xnode,xG)
x=calc_x(elem,xnode,xi);
I = (x-xG)'*(x-xG);
return
 
function n=quadratureorder(elem)
n=2*orderN(elem);
return