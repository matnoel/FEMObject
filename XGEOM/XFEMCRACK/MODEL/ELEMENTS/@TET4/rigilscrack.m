function ke = rigilscrack(elem,node,ls,varargin)

if ~isa(ls,'LSCRACK')
    error('l''argument doit etre une LSCRACK')
end

mat = getmaterial(elem);

ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
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

switch getlstype(elem)
    case {'indomain','in'}
        ls1 = lssupportval(connece);
        
        if all(ls1<=0) || all(ls1>=0)
        gauss=calc_gauss(elem,'rigi');       
        else
        [subgaussin,subgaussout] = calc_lssubgauss(eleme,ls1,'rigi');  
        gauss.coord = concat(subgaussin.coord,subgaussout.coord,4);
        gauss.w = concat(subgaussin.w,subgaussout.w,4);
        gauss.nbgauss = subgaussin.nbgauss + subgaussout.nbgauss;
        end
        
        phi = getN(eleme,gauss.coord)*ls1;
        Hphi = (phi>=0)-(phi<0);
            
         Bls = zerosND(3,getnbddl(eleme),1,gauss.nbgauss);   
         [B,detJ] = calc_B(eleme,xnodee,gauss.coord);
         for i=1:3
         repddli = 2*(i-1)+[1:2];
         Bi = B(:,repddli); 
             if connecenriche(i)
               Blsi = [Bi , Hphi*Bi];
             else
                 Blsi=Bi;
             end
         repddlienrich = sum(connecnbddle(1:i-1))+[1:connecnbddle(i)];
         Bls(:,repddlienrich)=Blsi;
         end
         D=calc_opmat(mat,eleme,xnodee,gauss.coord);
         ke(:,:,e) = sum(gauss.w*abs(detJ)*(Bls'*D*Bls),4);

    case 'cut'
    lse = lssupportval(connece);
    [subgaussin,subgaussout] = calc_lssubgauss(eleme,lse,'rigi');
    
    [B,detJ] = calc_B(eleme,xnodee,subgaussin.coord);
    Bls = assembleBls(B,-1);
    D = calc_opmat(mat,eleme,xnodee,subgaussin.coord);
    ke(:,:,e) = ke(:,:,e)+...
        sum(subgaussin.w*abs(detJ)*(Bls'*D*Bls),4);
    [B,detJ] = calc_B(eleme,xnodee,subgaussout.coord);
    Bls = assembleBls(B,1);
    D = calc_opmat(mat,eleme,xnodee,subgaussout.coord);
    ke(:,:,e) = ke(:,:,e)+...
        sum(subgaussout.w*abs(detJ)*(Bls'*D*Bls),4);
    case 'bicut'
        
        

    
end


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

