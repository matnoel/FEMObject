function fe = rigitangu(elem,node,q,U,varargin)

n = getcharin('intorder',varargin,'rigitang');
xnode = node(elem);
qe = localize(elem,q);
Ue = localize(elem,U);
fe = integrate(elem,xnode,n,@eval_feint,qe,Ue);

function fe = eval_feint(xi,elem,xnode,qe,Ue)
mat = getmaterial(elem);
fe = rigitangu(mat,elem,xnode,xi,qe,Ue);
return