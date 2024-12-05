function fe = bodysource(elem,node,ls,varargin)


switch getlstype(elem)
    case {'out'}
fe = zerosND(getnbddl(elem),1,getnbelem(elem));      
    otherwise
fe = zerosND(getnbddl(elem),1,getnbelem(elem));
xnode = getcoord(node,elem);
ls = getvalue(ls);
connec = getconnec(elem);

for e=1:getnbelem(elem)
    eleme = getelem(elem,e);
    connece = connec(e,:);
    lse = ls(connece);
     xnodee = xnode(:,:,e);    
    
    if all(lse<=0) 
    gauss = calc_gauss(eleme,2);
    elseif ~all(lse>=0)    
    [gauss,gaussout] = calc_lssubgauss(eleme,lse,2);  
    end
        N = calc_N(eleme,xnodee,gauss.coord);
        detJ = calc_detJ(eleme,xnodee,gauss.coord);
       fe(:,:,e)= sum(gauss.w*(N')*abs(detJ),4);
    
end
end

