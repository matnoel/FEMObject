function ke = masslsinterface_elem(elem,node,ls,varargin)

method = getcharin('method',varargin,'penal');

switch method
    case 'penal'
ke=zerosND(getnbddl(elem),getnbddl(elem));
xnodeelocal = nodelocalcoord(elem);
xnodeseg = contour_oneelem(ls,xnodeelocal);
elemseg = SEG2(xnodeseg,1,1:2);
    
gauss = calc_gauss(elemseg,2);
xi = calc_x(elemseg,xnodeseg,gauss.coord);
xnodeseg=double(xnodeseg);
x1 = calc_x(elem,node,xnodeseg(1,:));
x2 = calc_x(elem,node,xnodeseg(2,:));
    
le = norm(x2-x1);
N = calc_N(elem,node,xi);
    
ke(:,:)=sum(gauss.w*(le/2)*N'*N,4);


    case 'nitsche'
fact = getcharin('fact',varargin,1000);   
ke=zerosND(getnbddl(elem),getnbddl(elem));
xnodeelocal = nodelocalcoord(elem);
xnodeseg = contour_oneelem(ls,xnodeelocal);
elemseg = SEG2(xnodeseg,1,1:2);
gauss = calc_gauss(elemseg,2);
xi = calc_x(elemseg,xnodeseg,gauss.coord);
xnodeseg=double(xnodeseg);
x1 = calc_x(elem,node,xnodeseg(1,:));
x2 = calc_x(elem,node,xnodeseg(2,:));
h = rayon_du_cercle_inscrit(node);
le = norm(x2-x1);
N = calc_N(elem,node,xi);
DN = calc_DN(elem,node,xi);
    
ke(:,:)=-sum(gauss.w*(le/2)*DN'*(DN*ls)*N,4)'-sum(gauss.w*(le/2)*DN'*(DN*ls)*N,4)+(fact)*sum(gauss.w*(le/2)*N'*N,4); 

end