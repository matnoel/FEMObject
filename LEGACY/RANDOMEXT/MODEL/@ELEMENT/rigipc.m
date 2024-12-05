function ke = rigipc(elem,node,PC,varargin)

n=getcharin('intorder',varargin,'rigi');
xnode = node(elem);
gauss=calc_gauss(elem,n);
detJ=calc_detJ(elem,xnode,gauss.coord);
mat = getmaterial(elem);
ke = rigipc(mat,elem,xnode,gauss.coord,PC);
ke = sum(gauss.w*abs(detJ)*ke,4);

