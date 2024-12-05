function ke = rigilspcdomain(elem,node,PC,ls,varargin)
switch getlstype(elem)
    case {'in','indomain'}
     ke = rigipc(elem,node,PC,varargin{:});          
    case 'out'
    ke = zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
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
Nxk = Ntet4(gauss_spatial.coord);

mat = getmaterial(elem);
ke=zeros(getnbddl(elem),getnbddl(elem),getnbelem(elem),length(PC));

D = calc_opmat(mat,elem);
xnode = getcoord(node,elem);
[B,detJ] = calc_B(elem,xnode,[]);
Z = double(B'*D*B*detJ);

connec = getconnec(elem);
xnode = getcoord(node,elem);

lsval = getvalue(ls);
fprint('calcul elementaire')
for e=1:getnbelem(elem);
    pourcentage(e,getnbelem(elem),10)
    eleme= getelem(elem,e);
    connece = getconnec(eleme);
    xnodee = xnode(:,:,e);
    lse=double(lsval(connece));

    lsxk = Nxk*lse;
    Gmoins = (lsxk<0);
    %Halphaw = (diag(gauss_sto.w)*(Halpha));
    W=spdiags(gauss_sto.w(:),0,length(gauss_sto.w),length(gauss_sto.w));
    Halphaw = W*Halpha;
    %Halphaw = Halpha;
    %for ii=1:size(Halpha,1)
    %    Halphaw(ii,:)=gauss_sto.w(ii)*Halpha(ii,:);
    %end
    
    C = (gauss_spatial.w(:)'*Gmoins)*Halphaw;
    for alpha=1:length(PC)
    ke(:,:,e,alpha) = C(alpha)*Z(:,:,e);
    end    
end

ke = PCMYDOUBLEND(MYDOUBLEND(ke),PC,4);


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

function Bls = assembleBls(B,psi);

    Bls = [B ,factnodeB(B,psi)];   
    rep = zeros(3,4,2);
    rep(:)=1:numel(rep);
    rep = permute(rep,[1,3,2]);
    Bls(:,:) = Bls(:,rep(:)) ; 

return


function AKalpha = permute_ddl(AKalpha)
rep=[1,4,2,5,3,6];


AKalpha = AKalpha(rep,rep);

return
