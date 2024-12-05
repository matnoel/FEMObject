function ke = rigilsmaterial(elem,node,ls,varargin)

matin = getmaterial(ls);
matout = getmaterial(elem);

if ~isenrich(elem) %%% cas non enrichi     
switch getlstype(elem)
case {'indomain'}
ke = rigi(elem,node,varargin{:});
case 'in'
ke = rigi(setmaterial(elem,matin),node,varargin{:});
case 'cut'
ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
Din = calc_opmat(matin,elem);
Dout = calc_opmat(matout,elem);
connec = getconnec(elem);
xnode = getcoord(node,elem);
[B,detJ] = calc_B(elem,xnode,[]);
Zin = B'*Din*B*detJ;
Zout = B'*Dout*B*detJ;
lsval = getvalue(ls);
for e=1:getnbelem(elem)

    connece = connec(e,:);
    
    lse = lsval(connece);    
    if all(lse<=0) 
    ke(:,:,e)=Zin(:,:,e)*1/2;
    elseif all(lse>=0) 
    ke(:,:,e)=Zout(:,:,e)*1/2;    
    else
    [subgaussin,subgaussout] = elem_lssubgauss(elem,lse,0,1);
    ke(:,:,e)=Zin(:,:,e)*sum(subgaussin.w)+Zout(:,:,e)*sum(subgaussout.w);
    end
end
otherwise
        error('type non prevu')
end

else %%%% ca enrichi
switch getenrichtype(ls)
    case 3
n=getcharin('intorder',varargin,4);        
   otherwise
n=getcharin('intorder',varargin,2);
end

switch getlstype(elem)
    case {'in','indomain'}
        error('pas programme : cas d''un enrichissement partiel de l''element')
    case {'cut','touchcut'}
ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
lsval = getvalue(ls);

for e=1:getnbelem(elem)
    eleme = getelem(elem,e);
    connece = getconnec(eleme);
    xnodee = getcoord(node,getconnec(eleme)');
    lse = lsval(connece);
    
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
        error('type non prevu')

end
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


