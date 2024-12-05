function [fe,L] = fintpc(elem,node,PC,q,varargin)

n=getcharin('intorder',varargin,'rigitang');
xnode = node(elem);
qe=localizepc(elem,q);
xnode = node(elem);
gauss=calc_gauss(elem,n);
detJ=calc_detJ(elem,xnode,gauss.coord);
mat = getmaterial(elem);
fe = fintpc(mat,elem,xnode,gauss.coord,qe,PC);
fe = sum(gauss.w*abs(detJ)*fe,4);

