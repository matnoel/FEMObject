function ke = rigilscrack_tip(elem,node,ls,varargin)

tipnumber = getparam(elem,'tipnumber');

lstip = getlstip(ls,tipnumber); 
lstipval = getvalue(lstip);
lssupport=getlssupport(ls);
lssupportval = getvalue(lssupport);

orderquad=2;

mat = getmaterial(elem);

ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
connec = getconnec(elem);
connecenrich = getparam(elem,'connecenrich');
connecnbddl = getparam(elem,'connecnbddl');
conneclsenrichnature = getlsenrichnature(node,connec);
xnode = getcoord(node,elem);


for e=1:getnbelem(elem);
    eleme= getelem(elem,e);
    connece = getconnec(eleme);
    nodee = getnode(node,connece);
    connecenriche = connecenrich(e,:);
    connecnbddle = connecnbddl(e,:);
    conneclsenrichnaturee = conneclsenrichnature(e,:);
    xnodee = xnode(:,:,e);

        ls1 = lssupportval(connece);
        ls2 = lstipval(connece);
          
        if (all(ls1<=0) || all(ls1>=0)) && (all(ls2<=0) || all(ls2>=0)) 
        gauss=calc_gauss(elem,orderquad);       
        elseif (all(ls1<=0) || all(ls1>=0))
        [subgaussin,subgaussout] = calc_lssubgauss(eleme,ls2,orderquad);  
        gauss.coord = concat(subgaussin.coord,subgaussout.coord,4);
        gauss.w = concat(subgaussin.w,subgaussout.w,4);
        gauss.nbgauss = subgaussin.nbgauss + subgaussout.nbgauss;
        elseif (all(ls2<=0) || all(ls2>=0))
        [subgaussin,subgaussout] = calc_lssubgauss(eleme,ls1,orderquad);  
        gauss.coord = concat(subgaussin.coord,subgaussout.coord,4);
        gauss.w = concat(subgaussin.w,subgaussout.w,4);
        gauss.nbgauss = subgaussin.nbgauss + subgaussout.nbgauss;
        else      
        [subgaussinin,subgaussinout,subgaussoutin,subgaussoutout] = calc_bilssubgauss(eleme,ls1,ls2,orderquad);  
%        gaussorderspatial = getcharin('gaussorderspatial',varargin,2);
%        nbsubgaussspatial =    getcharin('nbsubgaussspatial',varargin,3) ;
%        gauss = calc_lssubgauss(elem,gaussorderspatial,nbsubgaussspatial);
        gauss.coord = subgaussinin.coord;
        gauss.coord = concat(gauss.coord,subgaussinout.coord,4);
        gauss.coord = concat(gauss.coord,subgaussoutin.coord,4);
        gauss.coord = concat(gauss.coord,subgaussoutout.coord,4);
        gauss.w = subgaussinin.w;
        gauss.w = concat(gauss.w,subgaussinout.w,4);
        gauss.w = concat(gauss.w,subgaussoutin.w,4);
        gauss.w = concat(gauss.w,subgaussoutout.w,4);
        gauss.nbgauss = size(gauss.w,4);
        end
        
         Bls = zerosND(3,getnbddl(eleme),1,gauss.nbgauss);   
         [B,detJ,DN] = calc_B(eleme,xnodee,gauss.coord);
         N = calc_N(eleme,xnodee,gauss.coord);
         
        phi = getN(eleme,gauss.coord)*ls1;
        Hphi = (phi>=0)-(phi<0);
        phi2 = getN(eleme,gauss.coord)*ls2;
        isphi2neg = phi2<=0;
        psi = -Hphi.*phi2.*isphi2neg;       
        dpsi = Hphi.*(-DN*ls2).*isphi2neg;
        Dpsi = zerosND([3,2,sizeND(dpsi)]);
        Dpsi(1,1) = dpsi(1);
        Dpsi(2,2) = dpsi(2);
        Dpsi(3,1) = dpsi(2);
        Dpsi(3,2) = dpsi(1);
   
         for i=1:3
         repddli = 2*(i-1)+[1:2];
         Bi = B(:,repddli); 
         Ni = N(:,repddli);
         
             if connecenriche(i) && strncmp(conneclsenrichnaturee{i},'sup',3)
               Blsi = [Bi , Hphi*Bi];
             elseif connecenriche(i) && strncmp(conneclsenrichnaturee{i},'tip',3)
               Blsi = [Bi , psi*Bi + Dpsi*Ni];  
             else
               Blsi=Bi;
             end
         repddlienrich = sum(connecnbddle(1:i-1))+[1:connecnbddle(i)];
         Bls(:,repddlienrich)=Blsi;
         end
         D=calc_opmat(mat,eleme,xnodee,gauss.coord);
         ke(:,:,e) = sum(gauss.w*abs(detJ)*(Bls'*D*Bls),4);
        
               
 
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

