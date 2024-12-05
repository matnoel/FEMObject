function ke = rigitang(elem,node,q,varargin)
% function ke = rigitang(elem,node,q,varargin)

n = getcharin('intorder',varargin,'rigitang');
xnode = node(elem);
qe = localize(elem,q);

ke = integrate(elem,xnode,n,@eval_ketang,qe);

function ke = eval_ketang(xi,elem,xnode,qe)
mat = getmaterial(elem);
ke = rigitang(mat,elem,xnode,xi,qe);
return