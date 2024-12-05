function ke = rigils(elem,node,ls,varargin)
ls=getlevelset(ls,1);
isin = lsisin(elem,ls);
isout = lsisout(elem,ls);
iscut = lsiscut(elem,ls);
elemin = getelem(elem,isin);
elemout = getelem(elem,isout);
elemcut = getelem(elem,iscut);

if isa(ls,'LEVELSETMAT') 
    matin = getmaterial(ls);
    matout = getmaterial(elem);
else
    matin = getmaterial(elem);
    matout = getmaterial(ls);
end

ke=zerosND(elem.nbddl,elem.nbddl,getnbelem(elem));

n=getcharin('intorder',varargin,'rigi');
gauss=calc_gauss(elem,n);

if getnbelem(elemin) & ~isempty(matin)
    xnodein = getcoord(node,getconnec(elemin)');
for l=1:gauss.nbgauss
D=calc_opmat(matin,elemin,xnodein,gauss.coord(l,:));
if getenrich(ls)==2
  [B,detJ]=calc_Bls(elemin,xnodein,gauss.coord(l,:),ls,'in');  
else
  [B,detJ]=calc_B(elemin,xnodein,gauss.coord(l,:));
end
ke(:,:,isin)=ke(:,:,isin)+gauss.w(l)*abs(detJ)*B'*D*B;
end
end

if getnbelem(elemout) & ~isempty(matout)
xnodeout = getcoord(node,getconnec(elemout)');
for l=1:gauss.nbgauss
D=calc_opmat(matout,elemout,xnodeout,gauss.coord(l,:));
if getenrich(ls)==2
  [B,detJ]=calc_Bls(elemout,xnodeout,gauss.coord(l,:),ls,'out');  
else
  [B,detJ]=calc_B(elemout,xnodeout,gauss.coord(l,:));
end
ke(:,:,isout)=ke(:,:,isout)+gauss.w(l)*abs(detJ)*B'*D*B;
end
end

if getnbelem(elemcut)>0
n=getcharin('intorder',varargin,'rigils');
gauss=calc_gauss(elem,n);   
xnodecut = getcoord(node,getconnec(elemcut)');
[xnodein,xnodeout,elemin,elemout,nodeplus]=lsdivideelem(elemcut,node,ls);

for l=1:gauss.nbgauss
if ~isempty(matin)
gaussin = calc_x(elemin,xnodein,gauss.coord(l,:));
detJin = calc_detJ(elemin,xnodein,gauss.coord(l,:));
D=calc_opmat(matin,elemcut,xnodecut,gaussin);
if getenrich(ls)>0
  [B,detJ]=calc_Bls(elem,xnodecut,gaussin,ls,'in');
else
  [B,detJ]=calc_B(elem,xnodecut,gaussin);
end
ke(:,:,iscut)=ke(:,:,iscut)+gauss.w(l)*abs(detJ)*abs(detJin)*B'*D*B;
end

if ~isempty(matout)
gaussout = calc_x(elemout,xnodeout,gauss.coord(l,:));
detJout = calc_detJ(elemout,xnodeout,gauss.coord(l,:));
D=calc_opmat(matout,elemcut,xnodecut,gaussout);
if getenrich(ls)>0
[B,detJ]=calc_Bls(elem,xnodecut,gaussout,ls,'out');
else
  [B,detJ]=calc_B(elem,xnodecut,gaussout);
end
ke(:,:,iscut)=ke(:,:,iscut)+gauss.w(l)*abs(detJ)*abs(detJout)*B'*D*B;
end
end

end

