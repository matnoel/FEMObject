function ke = rigils(elem,node,ls,varargin)
ls=getlevelset(ls,1);
isin = lsisin(elem,ls);
isout = lsisout(elem,ls);
iscut = lsiscut(elem,ls);

if isa(ls,'LEVELSETMAT') 
    matin = getmaterial(ls);
    matout = getmaterial(elem);
else
    matin = getmaterial(elem);
    matout = getmaterial(ls);
end

if all(isin) | all(isout) 
n=getcharin('intorder',varargin,'rigi');
gauss=calc_gauss(elem,n);   
xnode = getcoord(node,getconnec(elem)');

ke=zerosND(elem.nbddl,elem.nbddl,getnbelem(elem));

if all(isin)
mat=matin ;
choix = 'in';
else
mat=matout;
choix = 'out';
end

if ~isempty(mat)
    for l=1:gauss.nbgauss
D=calc_opmat(mat,elem,xnode,gauss.coord(l,:));
if getenrich(ls)==2
  [B,detJ]=calc_Bls(elem,xnode,gauss.coord(l,:),ls,choix);  
else
  [B,detJ]=calc_B(elem,xnode,gauss.coord(l,:));
end
ke=ke+gauss.w(l)*abs(detJ)*B'*D*B;
    end
end

elseif all(iscut)
n=getcharin('intorder',varargin,'rigils');
gauss=calc_gauss(elem,n);   
xnode = getcoord(node,getconnec(elem)');
ke=zerosND(elem.nbddl,elem.nbddl,getnbelem(elem));  

[xnodein,xnodeout,elemin,elemout,nodeplus]=lsdivideelem(elem,node,ls);


for l=1:gauss.nbgauss

gaussin = calc_x(elemin,xnodein,gauss.coord(l,:));
detJin = calc_detJ(elemin,xnodein,gauss.coord(l,:));
D=calc_opmat(matin,elem,xnode,gaussin);
if getenrich(ls)>0
    [B,detJ]=calc_Bls(elem,xnode,gaussin,ls,'in');
else
  [B,detJ]=calc_B(elem,xnode,gaussin);
end
ke=ke+gauss.w(l)*abs(detJ)*abs(detJin)*B'*D*B;


if ~isempty(matout)
gaussout = calc_x(elemout,xnodeout,gauss.coord(l,:));
detJout = calc_detJ(elemout,xnodeout,gauss.coord(l,:));
D=calc_opmat(matout,elem,xnode,gaussout);
if getenrich(ls)>0
[B,detJ]=calc_Bls(elem,xnode,gaussout,ls,'out');
else
  [B,detJ]=calc_B(elem,xnode,gaussout);
end
ke=ke+gauss.w(l)*abs(detJ)*abs(detJout)*B'*D*B;
end

end

else
    error('elements de types diffrents -> utiliser lssplitelem')
end