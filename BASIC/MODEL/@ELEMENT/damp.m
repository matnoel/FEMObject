function me = damp(elem,node,varargin)
% function me = damp(elem,node,varargin)

n=getcharin('intorder',varargin,'damp');
xnode = node(elem);

me = integrate(elem,xnode,n,@eval_me);


function me = eval_me(xi,elem,xnode)
mat = getmaterial(elem);
me = damp(mat,elem,xnode,xi);
return
