function ke = rigi(elem,node,varargin)
% function ke = rigi(elem,node,varargin)

mat = getmaterial(elem);
if isa(mat,'FOUR_ISOT') && isparam(mat,'r')
    n=getcharin('intorder',varargin,'mass');
else
    n=getcharin('intorder',varargin,'rigi');
end
xnode = node(elem);

ke = integrate(elem,xnode,n,@eval_ke,varargin{:});

function ke = eval_ke(xi,elem,xnode,varargin)
mat = getmaterial(elem);
ke = rigi(mat,elem,xnode,xi,varargin{:});
return