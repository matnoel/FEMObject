function ke = adv(elem,node,varargin)
% function ke = adv(elem,node,varargin)

n=getcharin('intorder',varargin,'adv');
xnode = node(elem);

ke = integrate(elem,xnode,n,@eval_ke,varargin{:});

function ke = eval_ke(xi,elem,xnode,varargin)
mat = getmaterial(elem);
ke = adv(mat,elem,xnode,xi,varargin{:});
return