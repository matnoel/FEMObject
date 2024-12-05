function ke = rigiexpect(elem,node,a,b,varargin)

n = getcharin('intorder',varargin,'rigi');
xnode = node(elem);
ke = integrate(elem,xnode,n,@eval_ke,a,b);

function ke = eval_ke(xi,elem,xnode,a,b)
mat=getmaterial(elem);
ke = rigiexpect(mat,elem,xnode,xi,a,b);
return