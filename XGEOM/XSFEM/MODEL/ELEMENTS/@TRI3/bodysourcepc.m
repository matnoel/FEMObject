
function fe = bodysourcepc(elem,node,PC,ls,varargin)

 switch getlstype(elem)
      case {'out'}
fe = zerosND(getnbddl(elem),1,getnbelem(elem));      
    otherwise
        
fe=zeros(getnbddl(elem),1,getnbelem(elem),length(PC));

nbgausssto = getcharin('nbgausssto',varargin,getorder(PC)+2);
nbsubgausssto =  getcharin('nbsubgausssto',varargin,4); %pour les fonc non regulieres.
type=getpoly(PC,1);
%
type_filtre = getcharin('type_filtre',varargin);
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
    
gauss_sto = calc_gausslobattopoints(type,nbgausssto);
end

ls = randomeval(ls,gauss_sto.coord,RANDVARS(PC));%calcul de ls dans les pts de gauss stoch.
Halpha = polyval(PC,gauss_sto.coord); 
w = MYDOUBLEND(reshape(gauss_sto.w,[1,1,length(gauss_sto.w)]));
Halpha = MYDOUBLEND(reshape(full(Halpha),[1,1,length(gauss_sto.w),length(PC)]));

xnode = getcoord(node,elem);

method = getcharin('method',varargin);
 switch method 
        case {'penal','pas de methode'}
for e=1:getnbelem(elem);
    pourcentage(e,getnbelem(elem),10)
    fetemp = zerosND(getnbddl(elem),1,gauss_sto.nbgauss);
    eleme= getelem(elem,e);
    connece = getconnec(eleme);
    xnodee = xnode(:,:,e);
    lse=ls(connece);
    for i=1:gauss_sto.nbgauss
    [subgaussin,subgaussout] = elem_lssubgauss(eleme,lse{i},6); 
    detJ = calc_detJ(eleme,xnodee,subgaussin.coord);
    Nxk = Ntri3(subgaussin.coord);
    fetemp(:,:,i) = double(subgaussin.w(:)'*((Nxk))*abs(detJ))';
    end
 
    fetemp = (fetemp*w)*Halpha;
    fe(:,:,e,:) = sum(fetemp,3);

end
fe = PCMYDOUBLEND(MYDOUBLEND(fe),PC,4);

     case 'nitsche'
for e=1:getnbelem(elem);
    pourcentage(e,getnbelem(elem),10)
    fetemp = zerosND(getnbddl(elem),1,gauss_sto.nbgauss);
    eleme= getelem(elem,e);
    connece = getconnec(eleme);
    xnodee = xnode(:,:,e);
    lse=ls(connece);
    for i=1:gauss_sto.nbgauss
    [subgaussin,subgaussout] = elem_lssubgauss(eleme,lse{i},0,2); 
    detJ = calc_detJ(eleme,xnodee,subgaussin.coord);
    Nxk = Ntri3(subgaussin.coord);
    fetemp(:,:,i) = double(subgaussin.w(:)'*((Nxk))*abs(detJ))';
    end
 
    fetemp = (fetemp*w)*Halpha;
    fe(:,:,e,:) = sum(fetemp,3);

end
fe = PCMYDOUBLEND(MYDOUBLEND(fe),PC,4);
        case 'carac'
for e=1:getnbelem(elem);
    pourcentage(e,getnbelem(elem),10)
    fetemp = zerosND(getnbddl(elem),1,gauss_sto.nbgauss);
    eleme= getelem(elem,e);
    connece = getconnec(eleme);
    xnodee = xnode(:,:,e);
    lse=ls(connece);
    for i=1:gauss_sto.nbgauss
    [subgaussin,subgaussout] = calc_lssubgauss(eleme,lse{i},4); 
    detJ = calc_detJ(eleme,xnodee,subgaussin.coord);
    Nxk = Ntri3(subgaussin.coord);
    psi=-Nxk*lse{i};
%psi = (Nxk*lse{i}).^2;
% psi =  exp(-Nxk*lse{i})-1;
 
    fetemp(:,:,i) = double(sum(subgaussin.w(:)'*psi*((Nxk))*abs(detJ),4)');
    end
    fetemp = (fetemp*w)*Halpha;
    fe(:,:,e,:) = sum(fetemp,3);

end

fe = PCMYDOUBLEND(MYDOUBLEND(fe),PC,4);
  
     
     case 'filtre'
    l=getcharin('l',varargin,0.05);
    a=getcharin('a',varargin,2);
    c=getcharin('c',varargin,5);   
    for e=1:getnbelem(elem);
    pourcentage(e,getnbelem(elem),10)
    fetemp = zerosND(getnbddl(elem),1,gauss_sto.nbgauss);
    eleme= getelem(elem,e);
    connece = getconnec(eleme);
    xnodee = xnode(:,:,e);
    lse=ls(connece);
    for i=1:gauss_sto.nbgauss
    %[subgaussin,subgaussout] = calc_lssubgauss(eleme,lse{i},4); 
    subgaussin = calc_gauss(eleme,4);
    detJ = calc_detJ(eleme,xnodee,subgaussin.coord);
    Nxk = Ntri3(subgaussin.coord);
    psi=-Nxk*lse{i};
    Hm=heaviside_filtered(double(psi),l,a,c,'type_filtre',type_filtre);
    fetemp(:,:,i) = double(sum(subgaussin.w(:)'*Hm*((Nxk))*abs(detJ),4)');
    end
    fetemp = (fetemp*w)*Halpha;
    fe(:,:,e,:) = sum(fetemp,3);

end

fe = PCMYDOUBLEND(MYDOUBLEND(fe),PC,4);
    
         
 end
end
end

