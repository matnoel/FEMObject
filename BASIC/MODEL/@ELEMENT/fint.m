function fe = fint(elem,node,q,varargin)
% function fe = fint(elem,node,q,varargin)

n=getcharin('intorder',varargin,'rigitang');
xnode = node(elem);
qe = localize(elem,q);

fe = integrate(elem,xnode,n,@eval_feint,qe);

function fe = eval_feint(xi,elem,xnode,qe)
mat = getmaterial(elem);
fe = fint(mat,elem,xnode,xi,qe);
return