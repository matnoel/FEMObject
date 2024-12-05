function ke = rigilsdomain(elem,node,ls,varargin)


switch getlstype(elem)
    case {'indomain','in'}
      ke = rigi(elem,node,varargin{:});
    case 'out'
      ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
    case 'cut'
     
ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
mat = getmaterial(elem);
D = calc_opmat(mat,elem);
xnode = getcoord(node,elem);

%if isaxi(elem)
%gausstemp = calc_gauss(elem,0);    
%[B,detJ] = calc_B(elem,xnode,gausstemp.coord);
%else
[B,detJ] = calc_B(elem,xnode,[]);    
%end

Z = B'*D*B*detJ;

ls = getvalue(ls);
connec = getconnec(elem);

for e=1:getnbelem(elem)
    connece = connec(e,:);
    lse = ls(connece);
    
    if all(lse<=0) 
    ke(:,:,e)=Z(:,:,e)*1/2;
    elseif ~all(lse>=0) 
    [subgaussin,subgaussout] = elem_lssubgauss(elem,lse,0,1);
    ke(:,:,e)=Z(:,:,e)*sum(subgaussin.w);
    
    end
end
        

    
    otherwise
        error('type non prevu')
        
        
end
        