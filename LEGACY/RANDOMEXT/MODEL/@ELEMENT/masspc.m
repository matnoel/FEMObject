function me = masspc(elem,node,PC,varargin)

n=getcharin('intorder',varargin,'mass');
xnode = node(elem);
gauss=calc_gauss(elem,n);
detJ=calc_detJ(elem,xnode,gauss.coord);
mat = getmaterial(elem);
me = masspc(mat,elem,xnode,gauss.coord,PC);
me = sum(gauss.w*abs(detJ)*me,4);
