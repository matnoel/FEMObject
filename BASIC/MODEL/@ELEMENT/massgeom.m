function me = massgeom(elem,node,varargin)
% function me = massgeom(elem,node,varargin)

n=getcharin('intorder',varargin,'mass');
xnode = getcoord(node,getconnec(elem)');

me = integrate(elem,xnode,n,@eval_me,varargin{:});

function me = eval_me(xi,elem,xnode,varargin)
N = calc_Ngeom(elem,xnode,xi,varargin{:});
me = N'*N;
return
