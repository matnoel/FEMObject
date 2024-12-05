function ke = rigilspcdomain(elem,node,PC,ls,varargin)
switch getlstype(elem)
    case {'in','indomain'}
     ke = rigipc(elem,node,PC,varargin{:});          
    case 'out'
    ke = zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
    case 'cut'
nbgausssto = getcharin('nbgausssto',varargin,getorder(PC)+1);
nbsubgausssto =  getcharin('nbsubgausssto',varargin,4);

type=getpoly(PC,1);
%
if  strcmp(class(type),'POLYFE')
gauss_sto = calc_gausspoints(PC,nbgausssto);
elseif  strcmp(class(type),'POLYLEGENDRE')
 if nbsubgausssto>1
 gauss_sto = calc_subgausspoints(PC,nbgausssto,nbsubgausssto);
 else
 gauss_sto = calc_gausspoints(PC,nbgausssto);
 end
elseif strcmp(class(type),'POLYFELAGRANGE')
    %gauss_sto.coord = getpoints(type);
    %gauss_sto.w = getweights(type);
    %gauss_sto.nbgauss = getnbpoints(type)
gauss_sto = calc_gausslobattopoints(type,nbgausssto);
end

gaussorderspatial = getcharin('gaussorderspatial',varargin,2);
subgausssto =    getcharin('nbsubgaussspatial',varargin,3) ;
gauss_spatial = elem_subgauss(elem,gaussorderspatial,subgausssto);

ls = randomeval(ls,gauss_sto.coord,RANDVARS(PC));
Halpha = polyval(PC,gauss_sto.coord); 
Nxk = Ntri3(gauss_spatial.coord);

mat = getmaterial(elem);
ke=zeros(getnbddl(elem),getnbddl(elem),getnbelem(elem),length(PC));

D = calc_opmat(mat,elem);
xnode = getcoord(node,elem);
%if isaxi(elem)
%gausstemp = calc_gauss(elem,0);    
%[B,detJ] = calc_B(elem,xnode,gausstemp.coord);
%else
[B,detJ] = calc_B(elem,xnode,[]);    
%end

Z = double(B'*D*B*detJ);

connec = getconnec(elem);
xnode = getcoord(node,elem);

lsval = getvalue(ls);

for e=1:getnbelem(elem);
    pourcentage(e,getnbelem(elem),10)
    eleme= getelem(elem,e);
    connece = getconnec(eleme);
    xnodee = xnode(:,:,e);
    lse=double(lsval(connece));

    lsxk = Nxk*lse;
    Gmoins = (lsxk<0);
    W=spdiags(gauss_sto.w(:),0,length(gauss_sto.w),length(gauss_sto.w));
    Halphaw = W*Halpha;
    
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
    rep = zeros(2,3,2);
    rep(:)=1:numel(rep);
    rep = permute(rep,[1,3,2]);
    Bls(:,:) = Bls(:,rep(:)) ; 

return


function AKalpha = permute_ddl(AKalpha)
rep=[1,4,2,5,3,6];


AKalpha = AKalpha(rep,rep);

return
