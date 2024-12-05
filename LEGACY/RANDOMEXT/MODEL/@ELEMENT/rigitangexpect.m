function ke = rigitangexpect(elem,node,q,a,b,varargin)

n = getcharin('intorder',varargin,'rigitang');
xnode = node(elem);
qe = localize(elem,q);
ke = integrate(elem,xnode,n,@eval_ketang,qe,a,b);

function ke = eval_ketang(xi,elem,xnode,qe,a,b)
mat = getmaterial(elem);
ke = rigitangexpect(mat,elem,xnode,xi,qe,a,b);
return