function ke = reac(elem,node,varargin)
% function ke = reac(elem,node,varargin)

n=getcharin('intorder',varargin,'mass');
xnode = node(elem);

ke = integrate(elem,xnode,n,@eval_ke,varargin{:});

function ke = eval_ke(xi,elem,xnode,varargin)
mat = getmaterial(elem);
ke = reac(mat,elem,xnode,xi,varargin{:});
return