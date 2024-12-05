function ke = rigilspcmaterial(elem,node,PC,ls,varargin)
switch getlstype(elem)
    case {'indomain','out'}
     ke = rigipc(elem,node,PC,varargin{:}); 
    case 'in'
        elem=setmaterial(elem,getmaterial(ls));
        ke = rigipc(elem,node,PC,varargin{:}); 
    case 'cut'
nbgausssto = getcharin('nbgausssto',varargin,getorder(PC)+1);
nbsubgausssto =  getcharin('nbsubgausssto',varargin,4) ;
if nbsubgausssto>1
gauss_sto = calc_subgausspoints(PC,nbgausssto,nbsubgausssto);
else
gauss_sto = calc_gausspoints(PC,nbgausssto);
end
gaussorderspatial = getcharin('gaussorderspatial',varargin,2);
subgausssto =    getcharin('nbsubgaussspatial',varargin,3) ;
gauss_spatial = elem_subgauss(elem,gaussorderspatial,subgausssto);

ls = randomeval(ls,gauss_sto.coord,RANDVARS(PC));
Halpha = polyval(PC,gauss_sto.coord); 
Nxk = Ntri3(gauss_spatial.coord);

matout = getmaterial(elem);
matin = getmaterial(ls);

ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem),length(PC));

Din = calc_opmat(matin,elem);
Dout = calc_opmat(matout,elem);
xnode = getcoord(node,elem);
[B,detJ] = calc_B(elem,xnode,[]);
Zin = B'*Din*B*abs(detJ);
Zout = B'*Dout*B*abs(detJ);

connec = getconnec(elem);
xnode = getcoord(node,elem);

lsval = getvalue(ls);

for e=1:getnbelem(elem);
    eleme= getelem(elem,e);
    connece = getconnec(eleme);
    xnodee = xnode(:,:,e);
    lse=double(lsval(connece));

    lsxk = Nxk*lse;
    Gin = (lsxk<0);
    Gout = lsxk>=0;
    Cin = (gauss_spatial.w(:)'*Gin)*(diag(gauss_sto.w)*(Halpha));
    Cout = (gauss_spatial.w(:)'*Gout)*(diag(gauss_sto.w)*(Halpha));
    for alpha=1:length(PC)
    ke(:,:,e,alpha) = Cin(alpha)*Zin(:,:,e)+Cout(alpha)*Zout(:,:,e);
    end    
end

ke = PCMYDOUBLEND(ke,PC,4);


end


function B = factnodeB(B,f)
if numel(f)==1
    B = B*f;
else
f = repmat(f(:)',2,1);
f = repmat(f(:)',size(B,1),1);
B = f.*B;
end
return

function Bls = assembleBls(B,psi)

    Bls = [B ,factnodeB(B,psi)];   
    rep = zeros(2,3,2);
    rep(:)=1:numel(rep);
    rep = permute(rep,[1,3,2]);
    Bls(:,:) = Bls(:,rep(:)) ; 

return


function AKalpha = permute_ddl(AKalpha)
rep=[1,4,2,5,3,6];


AKalpha = AKalpha(rep,rep);

return
