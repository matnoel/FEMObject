function ke = dispt(elem,node,varargin)
% function ke = dispt(elem,node,varargin)

n=getcharin('intorder',varargin,'rigi');
xnode = node(elem);

ke = integrate(elem,xnode,n,@eval_ke);

function ke = eval_ke(xi,elem,xnode)
mat = getmaterial(elem);
ke = dispt(mat,elem,xnode,xi);
return