function ke= rigilspcmaterial(elem,node,PC,ls,S,varargin)
% function ke= rigilspcmaterial(elem,node,PC,ls,S,PCls,varargin)
% --> calcul des matrices de rigidité élémentaire du groupe elem
% elem : groupe d'éléments finis
% node : noeuds correspondants
% PC : chaos utilisé pour la décomposition
% ls : level-set du model
% S : model
% PCls : PC utilisé pour la décomposition de ls
ls=getlevelset(ls,getlsnumber(elem));
%PCM = PCls;
%keyboard
%lspc = project(ls,PCM);
%lspc= lseval(lspc,S,PCM);
%lspcval = getvalue(lspc);
lspcval = getvalue(ls);
if ~isempty(getmaterial(ls))
    matin = getmaterial(ls);
    matout = getmaterial(elem);
else
    matin = getmaterial(elem);
    matout = getmaterial(ls);
end


n=getcharin('intorder',varargin,'rigils');
ng=getcharin('intordersto',varargin,getorder(PC)+1);

tic
switch getlstype(elem)
    case {'in','indomain','out'}
switch getlstype(elem)
    case {'in'}
        elem = setmaterial(elem,matin);
    case {'out'}
        elem = setmaterial(elem,matout);
end
if ~isempty(getmaterial(elem))
ke = rigi(elem,node);
ke = PCMYDOUBLEND(ke,one(PC),[]);
else
ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
end

    case 'cut'
       
ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem),length(PC));

for e=1:getnbelem(elem);
    fprintf('element %d/%d',e,getnbelem(elem));
    eleme = getelem(elem,e);
    xnode = getcoord(node,getconnec(eleme)');
    stomesh = getparam(eleme,'stomesh');
    gauss = calc_gausspoints(stomesh{1},ng);
    w = MYDOUBLEND(reshape(gauss.w,[1,1,1,length(gauss.w)]));
    xg = transfer(RANDVARS(stomesh{1}),RANDVARS(PC),gauss.coord);
    n = ceil(size(xg,1)/10);
    nb = floor(size(xg,1)/n); 
    H = [];
    for i=1:nb
       H = [H ; polyval(PC,xg((i-1)*n+1:i*n,:))]; 
    end   
    if size(xg,1)>n*nb
       H = [H ; polyval(PC,xg(n*nb+1:end,:))];
    end
    H = MYDOUBLEND(reshape(full(H),[1,1,1,length(gauss.w),length(PC)]));
    Bls = calc_B(eleme,xnode,gauss.coord(1,:));
    detJ = calc_detJ(eleme,xnode,gauss.coord(1,:));
    ketemp=zerosND(getnbddl(elem),getnbddl(elem),1,gauss.nbgauss);
    
    
    Din=calc_opmat(matin,eleme,xnode,gauss.coord(1,:));
    Ain = Bls'*Din*Bls*abs(detJ);
    if isempty(matout)
        Aout = 0;
    else
    Dout=calc_opmat(matout,eleme,xnode,gauss.coord(1,:));
    Aout = Bls'*Dout*Bls*abs(detJ);
    end
    rep = find(gauss.state==1);
    if ~isempty(rep)
    if ~isempty(matout)
    ketemp(:,:,1,rep) = Aout*(1/2);
    end
    end
    rep = find(gauss.state==-1);
    if ~isempty(rep)
    ketemp(:,:,1,rep) = Ain*(1/2);
    end
    rep = find(gauss.state==0);
    connec1 = getconnec(eleme);
    lspce = lspcval(connec1);
    lsk = randomlimiteval(lspce,gauss.coord(rep,:),RANDVARS(stomesh));
    lsk = MULTIMATRIX(lsk,size(lspce),[length(rep),1]);
    for j =1:length(rep)
    [subgaussin,subgaussout] = lssubgauss_oneelem(lsk{j},0,1);
    ketemp(:,:,1,rep(j)) = Ain*sum(subgaussin.w)+Aout*sum(subgaussout.w);
    end
if size(ketemp,4)<6000 && length(PC)<30
ketemp = (ketemp*w)*H;
ketemp = sum(ketemp,4);
elseif size(ketemp,4)>6000 && length(PC)<30
a = ceil(size(ketemp,4)/6000);
s = ceil(size(ketemp,4)/a);
ketemp1 = 0;
for p=1:a-1
   ketemp1 = ketemp1 + sum(ketemp(:,:,1,((p-1)*s)+1:p*s)*w(:,:,:,((p-1)*s)+1:p*s)*H(:,:,:,((p-1)*s)+1:p*s,:),4);
end
ketemp1 = ketemp1 + sum(ketemp(:,:,1,((a-1)*s)+1:end)*w(:,:,:,((a-1)*s)+1:end)*H(:,:,:,((a-1)*s)+1:end,:),4);
ketemp = ketemp1;
elseif length(PC)>30
    s = ceil(size(ketemp,4)/8);
    a = ceil(size(ketemp,4)/s);
    ketemp1 = 0;
for p=1:a-1
   ketemp1 = ketemp1 + sum(ketemp(:,:,1,((p-1)*s)+1:p*s)*w(:,:,:,((p-1)*s)+1:p*s)*H(:,:,:,((p-1)*s)+1:p*s,:),4);
end
ketemp1 = ketemp1 + sum(ketemp(:,:,1,((a-1)*s)+1:end)*w(:,:,:,((a-1)*s)+1:end)*H(:,:,:,((a-1)*s)+1:end,:),4);
ketemp = ketemp1;
end
ke(:,:,e,:) = permute(ketemp,[1,2,3,5,4]);        
end

ke = PCMYDOUBLEND(ke,PC,4);

end

fprintf('Elapsed time is %.3f seconds.',toc)