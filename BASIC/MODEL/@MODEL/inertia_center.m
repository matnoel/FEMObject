function s=inertia_center(S,varargin)

s=integrate(S,@quadorder,@fun);
s=s/measure(S);

function I = fun(xi,elem,xnode)
I=calc_x(elem,xnode,xi);
return

function n=quadorder(elem)
n=orderN(elem);
return