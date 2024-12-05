function ke = rigilspccrack(elem,node,PC,ls,varargin)

if ~isa(ls,'LSCRACK')
    error('l''argument doit etre une LSCRACK')
end


nbgausssto = getcharin('nbgausssto',varargin,getorder(PC)+1);
nbsubgausssto =  getcharin('nbsubgausssto',varargin,4) ;
if nbsubgausssto>1
gauss_sto = calc_subgausspoints(PC,nbgausssto,nbsubgausssto);
else
gauss_sto = calc_gausspoints(PC,nbgausssto);
end
gaussorderspatial = getcharin('gaussorderspatial',varargin,2);
nbsubgaussspatial =    getcharin('nbsubgaussspatial',varargin,3) ;
gauss_spatial = elem_subgauss(elem,gaussorderspatial,nbsubgaussspatial);


%Nsimul = 5000;
%gauss_sto.coord=2*rand(Nsimul,1)-1;
%gauss_sto.w=repmat(1/Nsimul,1,Nsimul);
%gauss_sto.nbgauss=Nsimul;

%gauss_sto = calc_gausspoints(PC,20);
%gauss_spatial = elem_subgauss(elem,2,3);

%w_spatial_sto = gauss_spatial.w(:)*gauss_sto.w(:)';
ls = randomeval(ls,gauss_sto.coord,RANDVARS(PC));
Halpha = polyval(PC,gauss_sto.coord); 
Nxk = Ntri3(gauss_spatial.coord);

mat = getmaterial(elem);
ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem),length(PC));

connec = getconnec(elem);
connecenrich = getparam(elem,'connecenrich');
connecnbddl = getparam(elem,'connecnbddl');
conneclsenrichnature = getlsenrichnature(node,connec);
xnode = getcoord(node,elem);

lssupport=getlssupport(ls);
lssupportval = getvalue(lssupport);


for e=1:getnbelem(elem);
    eleme= getelem(elem,e);
    connece = getconnec(eleme);
    nodee = getnode(node,connece);
    connecenriche = connecenrich(e,:);
    connecnbddle = connecnbddl(e,:);
    conneclsenrichnaturee = conneclsenrichnature(e,:);
    xnodee = xnode(:,:,e);

repddl=[];
for i=1:3
repi=2*(i-1)+[1:2];
repddl=[repddl,repi];
if connecenriche(i);
    repddl=[repddl,repi+6];
end
end
   
    
    switch getlstype(elem)
    
        case {'indomain','in'}
    [B,detJ] = calc_B(eleme,xnodee,[]);   
    D = calc_opmat(mat,eleme,xnodee,[]);
    lse = double(lssupportval(connece));
    lsxk = Nxk*lse;
    Gplus = (lsxk>=0);
    Gmoins = ~Gplus;
    C1 = 1/2*(gauss_sto.w*Halpha);
    C2 = (gauss_spatial.w(:)'*Gplus-gauss_spatial.w(:)'*Gmoins)*(diag(gauss_sto.w)*(Halpha));
    Z0 = detJ*B'*D*B;
    for alpha=1:length(PC)
    temp = [C1(alpha)*Z0,C2(alpha)*Z0;C2(alpha)*Z0,C1(alpha)*Z0];    
    temp = temp(repddl,repddl);
    ke(:,:,e,alpha) = temp;
    end    
  
         
        case 'cut'
    [B,detJ] = calc_B(eleme,xnodee,[]);   
    D = calc_opmat(mat,eleme,xnodee,[]);
    lse = double(lssupportval(connece));
    lsxk = Nxk*lse;
    Gplus = (lsxk>=0);
    Gmoins = ~Gplus;
    C1 = 1/2*(gauss_sto.w*Halpha);
    
    C2 = (gauss_spatial.w(:)'*Gplus-gauss_spatial.w(:)'*Gmoins)*(diag(gauss_sto.w)*(Halpha));
    Z0 = detJ*B'*D*B;
    
    for alpha=1:length(PC)
    temp = [C1(alpha)*Z0,C2(alpha)*Z0;C2(alpha)*Z0,C1(alpha)*Z0];    
    temp = temp(repddl,repddl);
    ke(:,:,e,alpha) = temp;
    end



        case 'bicut'
        
        
        
        otherwise
            error('type d''element pas prevu');
    end


end

ke = PCMYDOUBLEND(ke,PC,4);



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
