function ke = rigitangpc(elem,node,PC,q,varargin)

n=getcharin('intorder',varargin,'rigitang');
xnode = node(elem);
qe=localizepc(elem,q);

gauss=calc_gauss(elem,n);
detJ=calc_detJ(elem,xnode,gauss.coord);
mat = getmaterial(elem);
ke = rigitangpc(mat,elem,xnode,gauss.coord,qe,PC);
ke = sum(gauss.w*abs(detJ)*ke,4);

