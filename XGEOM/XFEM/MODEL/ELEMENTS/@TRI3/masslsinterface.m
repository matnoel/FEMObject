function ke = masslsinterface(elem,node,ls,varargin)


switch getlstype(elem)
    case {'in','indomain','out'}
ke = zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));      
    case 'cut'
ls = getlevelset(ls,getlsnumber(elem));   
ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
xnode = getcoord(node,elem);
ls = getvalue(ls);
connec = getconnec(elem);

for e=1:getnbelem(elem)
    eleme = getelem(elem,e);
    connece = connec(e,:);
    lse = ls(connece);
    xnodee = xnode(:,:,e);
    xnodeelocal = nodelocalcoord(elem);
    
    xnodeseg = contour_oneelem(lse,xnodeelocal);
    elemseg = SEG2(xnodeseg,1,1:2);
    gauss = calc_gauss(elemseg,2);
    xi = calc_x(elemseg,xnodeseg,gauss.coord);
    xnodeseg=double(xnodeseg);
    le = norm(xnodeseg(2,:)-xnodeseg(1,:));
    N = calc_N(eleme,xnodee,xi);
    ke(:,:,e)=sum(gauss.w*le*N'*N,4);
    
    
end
end

