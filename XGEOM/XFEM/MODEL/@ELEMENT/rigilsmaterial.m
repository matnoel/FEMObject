function ke = rigilsmaterial(elem,node,ls,varargin)

n=getcharin('intorder',varargin,'rigils');

matout = getmaterial(elem);
matin = getmaterial(ls);
switch getlstype(elem)
    case 'indomain'
      ke = rigi(elem,node,varargin{:});
    case 'in'
      ke = rigi(setmaterial(elem,matin),node,varargin{:});  
    case 'cut'
ke=zerosND(elem.nbddl,elem.nbddl,getnbelem(elem));
xnode = getcoord(node,elem);
for e=1:getnbelem(elem)
    eleme = getelem(elem,e);
    connece = getconnec(eleme);
    xnodee = xnode(:,:,e);
    lse = getvalue(ls);lse = lse(connece);
    
    if lsisout(eleme,lse)
        gauss = calc_gauss(eleme,n);
        ke(:,:,e) = ke(:,:,e) + integrate_with_gauss(eleme,xnodee,gauss,@eval_ke,matout,ls,'out');
    elseif lsisin(eleme,lse)
        gauss = calc_gauss(eleme,n);
        ke(:,:,e) = ke(:,:,e) + integrate_with_gauss(eleme,xnodee,gauss,@eval_ke,matin,ls,'in');
    elseif lsiscut(eleme,lse)
        [gaussin,gaussout] = calc_lssubgauss(eleme,lse,n);
        ke(:,:,e) = ke(:,:,e) + integrate_with_gauss(eleme,xnodee,gaussin,@eval_ke,matin,ls,'in');
        ke(:,:,e) = ke(:,:,e) + integrate_with_gauss(eleme,xnodee,gaussout,@eval_ke,matout,ls,'out');
    end
    
end
    otherwise
        error('pas prevu')
end


function ke = eval_ke(xi,elem,xnode,mat,ls,choix)

if getnbelem(elem)>0 
D=calc_opmat(mat,elem,xnode,xi);
Bls=calc_Bls(elem,xnode,xi,ls,choix);
ke = Bls'*D*Bls ;
else
ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
end


return

