function me = mass(elem,node,varargin)
% function me = mass(elem,node,varargin)

n=getcharin('intorder',varargin,'mass');
xnode = node(elem);

me = integrate(elem,xnode,n,@eval_me);

function me = eval_me(xi,elem,xnode)
mat = getmaterial(elem);
me = mass(mat,elem,xnode,xi);
return
