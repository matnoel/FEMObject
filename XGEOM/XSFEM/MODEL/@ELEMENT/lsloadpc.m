function fe= lsloadpc(elem,node,PC,S,compo,fun,varargin)
% compo : nom des composantes de la force
% fun : pointeur sur une fonction donnant la valeur des composantes en
% colonnes, peut N colonnes correspondant N cas de chargement
% appel de fun(x,varargin{:}) ou x sont le tableau de coordonnees de P points ou on veut evaluer la fonction
% size(x,1)=P size(x,2)=dimension de l'espace de definition des points
% la sortie est un 3D array de taille length(compo)*N*P

PC = getPC(PC);
if isa(S,'MODEL')
ls=getlevelset(S.ls,getlsnumber(elem));
else
ls=getlevelset(S,getlsnumber(elem));    
end
lspcval = getvalue(ls);
if ~isempty(getmaterial(ls))
    matin = getmaterial(ls);
    matout = getmaterial(elem);
else
    matin = getmaterial(elem);
    matout = getmaterial(ls);
end

n=getcharin('intorder',varargin,2*orderN(elem)); 
ng=getcharin('intordersto',varargin,getorder(PC)+1);

varargin = delcharin('intordersto',varargin);
xnode = getcoord(node,getconnec(elem)');

ftest=fun(double(xnode(1,:,1)),varargin{:});


switch getlstype(elem)
    case {'in','indomain','out'}
switch getlstype(elem)
    case {'in'}
        elem = setmaterial(elem,matin);
    case {'out'}
        elem = setmaterial(elem,matout);
end


if ~isempty(getmaterial(elem))
    fe = load(elem,node,compo,fun,varargin{:});
    %fe = repmat(fe,[1,1,1,numel(ls)]);
    fe = PCMYDOUBLEND(fe,one(PC),[]);
else

    fe = zerosND(elem.nbddl,size(ftest,2),getnbelem(elem));
end


    case 'cut'

    fe=zerosND(elem.nbddl,size(ftest,2),getnbelem(elem),length(PC));

    for e=1:getnbelem(elem)
        eleme = getelem(elem,e);
        xnode = getcoord(node,getconnec(eleme)');
        stomesh = getparam(eleme,'stomesh');
        gauss = calc_gausspoints(stomesh{1},ng);
        fetemp = zerosND(elem.nbddl,size(ftest,2),1,gauss.nbgauss);
        
        rep = find(gauss.state==1);
        if ~isempty(rep)
            eleme = setmaterial(eleme,matout);
        if ~isempty(matout)
        fef = load(eleme,node,compo,fun,varargin{:});
        fetemp(:,:,1,rep) = fef;
        end
        end
        
        rep = find(gauss.state==-1);
        if ~isempty(rep)
            eleme = setmaterial(eleme,matin);
            fef = load(eleme,node,compo,fun,varargin{:});
            fetemp(:,:,1,rep) = fef;
        end
        
        rep = find(gauss.state==0);
        connec1 = getconnec(eleme);
        lspce = lspcval(getpos(node,connec1));
        lsi = randomeval(lspce,gauss.coord(rep,:),RANDVARS(stomesh));
        lse = double(lsi);
        
        for j = 1:length(rep)
           I = sign(double(lse(:,j)));
        if all((I==1 | I==0))
           eleme = setmaterial(eleme,matout);
           if ~isempty(matout)
           fef = load(eleme,node,compo,fun,varargin{:});
           fetemp(:,:,1,rep(j)) = fef;
           end
        elseif all((I==-1 | I==0))
           eleme = setmaterial(eleme,matin);
           fef = load(eleme,node,compo,fun,varargin{:});
           fetemp(:,:,1,rep(j)) = fef;
        else 
       [subgaussin,subgaussout] = lssubgauss_oneelem(eleme,lsi{j},n,2);
       xiin = subgaussin.coord;
       xiin = reshape(xiin,[1 1 size(xiin,1),1]);
       xiout = subgaussout.coord;
       xiout = reshape(xiout,[1 1 size(xiout,1),1]);
       win = reshape(subgaussin.w,[1,1,size(subgaussin.w,1)]);
       wout = reshape(subgaussout.w,[1,1,size(subgaussout.w,1)]);
       detJ_in = calc_detJ(eleme,xnode,xiin);
       detJ_out = calc_detJ(eleme,xnode,xiout);
       fetemp(:,:,1,rep(j)) = fetemp(:,:,1,rep(j))  ...
           + sum(eval_fe(xiin,eleme,xnode,matin,compo,lsi{j},'in',...
              fun,varargin{:})*win*detJ_in,length(size(win)))...
           + sum(eval_fe(xiout,eleme,xnode,matout,compo,lsi{j},'out',...
              fun,varargin{:})*wout*detJ_out,length(size(wout)));
        end
        end

w = MYDOUBLEND(reshape(gauss.w,[1,1,1,length(gauss.w)]));

xg = transfer(RANDVARS(stomesh{1}),RANDVARS(PC),gauss.coord);
H = polyval(PC,xg);

H = MYDOUBLEND(reshape(full(H),[1,1,1,length(gauss.w),length(PC)]));
fetemp = (fetemp*w)*H;
fetemp = sum(fetemp,4);
fe(:,:,e,:) = permute(fetemp,[1,2,3,5,4]);


    end
    
    fe = PCMYDOUBLEND(fe,PC,4);
    
end
    

function fe = eval_fe(xi,elem,xnode,mat,compo,ls,choix,fun,varargin)
if ~isempty(mat) & getnbelem(elem)>0
[repddl,repddlcompo] = findddl(getddlnodedual(elem),compo);
N=calc_Nls(elem,xnode,xi,ls,choix);
x=calc_x(elem,xnode,xi);
f=fun(x,varargin{:});
f=MYDOUBLEND(f);
fe = N(repddl,:)'*f(repddlcompo,:);
else
ftest=fun(double(xnode(1,:)),varargin{:}) ;   
fe=zerosND(elem.nbddl,size(ftest,2),getnbelem(elem));
end

return

