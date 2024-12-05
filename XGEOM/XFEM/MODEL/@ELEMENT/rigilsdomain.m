function ke = rigilsdomain(elem,node,ls,varargin)

n=getcharin('intorder',varargin,'rigi');

switch getlstype(elem)
    case 'out'
      ke = zerosND(elem.nbddl,elem.nbddl,getnbelem(elem));      
    case {'in','indomain'}
      ke = rigi(elem,node,varargin{:});
    case 'cut'
ke=zerosND(elem.nbddl,elem.nbddl,getnbelem(elem));

for e=1:getnbelem(elem)
    eleme = getelem(elem,e);
    connece = getconnec(eleme);
    xnodee = getcoord(node,getconnec(eleme)');
    lse = getvalue(ls);lse = lse(connece);
    
    if lsisin(eleme,lse)
        gauss = calc_gauss(eleme,n);
        ke(:,:,e) = ke(:,:,e) + integrate_with_gauss(eleme,xnodee,gauss,@eval_ke);
    elseif lsiscut(eleme,lse)
        [gaussin,gaussout] = calc_lssubgauss(eleme,lse,n);
        ke(:,:,e) = ke(:,:,e) + integrate_with_gauss(eleme,xnodee,gaussin,@eval_ke);
    end
    
end
end


function ke = eval_ke(xi,elem,xnode)

if getnbelem(elem)>0 
D=calc_opmat(getmaterial(elem),elem,xnode,xi);
B=calc_B(elem,xnode,xi);
ke = B'*D*B ;
else
ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
end


return

