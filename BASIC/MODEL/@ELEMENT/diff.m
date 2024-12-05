function ke = diff(elem,node,varargin)
% function ke = diff(elem,node,varargin)

n=getcharin('intorder',varargin,'rigi');
xnode = node(elem);

ke = integrate(elem,xnode,n,@eval_ke,varargin{:});

function ke = eval_ke(xi,elem,xnode,varargin)
mat = getmaterial(elem);
ke = diff(mat,elem,xnode,xi,varargin{:});
return