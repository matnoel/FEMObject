function [conneccontour,nodecontour]=contour(elem,node,nodeval,contourval,varargin)


xnode=getcoord(node,getconnec(elem)');
dim=size(xnode,2);
[a,connec] = ismember(getconnec(elem),getnumber(node)) ;

nbelem = getnbelem(elem);

dim=size(xnode,2);
lsx=nodeval(connec)-contourval;
lsx=reshape(lsx,size(connec));
repzero=find(abs(lsx)<1e-10);
lsx(repzero)=0;


repcut = find(prod(lsx,2)<=0); % affichage meme ca coupe un noeud
elemcut = getelem(elem,repcut);
xnodecut = xnode(:,:,repcut);
xi = (lsx(repcut,1)+lsx(repcut,2))./ (lsx(repcut,1)-lsx(repcut,2)); 

xi=MYDOUBLEND(reshape(xi,1,1,length(repcut)));
x = calc_x(elemcut,xnodecut,xi);
nodecontour = double(x);
nodecontour = nodecontour(:);

conneccontour = [1:size(nodecontour,1)]';
nodecontour = NODE(nodecontour,conneccontour);

