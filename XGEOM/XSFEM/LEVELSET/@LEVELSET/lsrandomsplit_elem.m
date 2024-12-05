function [Ds,I_elemtot,stomesh] = lsrandomsplit_elem(ls,tol,anisotrope,statepoint,varargin)
% function [Ds,I_elemtot,stomesh] = lsrandomsplit_elem(ls,tol,PCls,anisotrope,statepoint,varargin)
% --> creation des partitions stochastiques de chaque element fini
% ls : LEVELSETS decomposee sur son chaos
% tol : tolerance k de decoupage
% PCls : chaos sur lequel la level-set sera decomposee
% anisotrope : 0 ou 1 respectivement decoupage isotrope ou anisotrope (pour
% anisotrope+ donner 2)
% statepoint : 0 ou 1 stockage de l'etat des sommets de chaque cellule dans
% le stomesh

D=getclassin('MODEL',varargin);
X = RANDVARS(ls);
stodim = getM(X);
if ~isa(getvalue(ls),'PCMATRIX')
    error('la level-set du model doit etre decomposee et evaluee sur son chaos')
end
lspcval = getvalue(ls);


U = RANDVARS();
for i=1:stodim
    U{i}=RVUNIFORM(0,1);
end
eoutglob = [];
einglob  = [];
ecutglob = [];

for p=1:getnbgroupelem(D)
    connec = getconnec(D.groupelem{p});
    numelem = getnumber(D.groupelem{p});
    fprintf('Group element %d ',p)
    
    for e=1:size(connec,1)
        fprintf('%d/%d ',e,size(connec,1))
        nume = numelem(e);
        cone = connec(e,:);
        lspce = lspcval(cone);
        
        if anisotrope==0
            I_elem{e}.order = zeros(0,stodim);
            I_elem{e}.way = zeros(0,stodim);
            
            xi=calc_sommets(calc_xi_ani(I_elem{e},stodim),stodim);
            I_elem{e}.signsommet = find_Iin(xi,stodim,U,X,lspce);
            
            I_elemtot{e} = cell(0,1);
            I_elemtot{e} = recursivecut_ND(I_elemtot{e},I_elem{e},tol,U,X,lspce,'stodim',stodim);
            
            stomesh{e} = POLYFEND(I_elemtot{e},'statepoint',statepoint);
        elseif anisotrope==1
            
            I_elem{e}.order = zeros(0,stodim);
            I_elem{e}.way = zeros(0,stodim);
            
            xi=calc_sommets(calc_xi_ani(I_elem{e},stodim),stodim);
            I_elem{e}.signsommet = find_Iin(xi,stodim,U,X,lspce);
            
            I_elemtot{e} = cell(0,1);
            I_elemtot{e} = recursivecut_ani_ND(I_elemtot{e},I_elem{e},tol,U,X,lspce,'stodim',stodim);
            
            stomesh{e} = POLYFEND(I_elemtot{e},'statepoint',statepoint);
        elseif anisotrope==2
            I_elem{e}.order = zeros(0,stodim);
            I_elem{e}.way = zeros(0,stodim);
            
            xi=calc_sommets(calc_xi_ani(I_elem{e},stodim),stodim);
            I_elem{e}.signsommet = find_Iin(xi,stodim,U,X,lspce);
            
            I_elemtot{e} = cell(0,1);
            I_elemtot{e} = recursivecut_ani_ND(I_elemtot{e},I_elem{e},tol,U,X,lspce,'stodim',stodim,'ani+');
            
            stomesh{e} = POLYFEND(I_elemtot{e},'statepoint',statepoint);
        end
        
        if I_elemtotin(I_elemtot{e})
            einglob=[einglob,nume];
        elseif I_elemtotout(I_elemtot{e})
            eoutglob=[eoutglob,nume];
        else
            ecutglob=[ecutglob,nume];
            
        end
        
    end
    
    D.groupelem{p} = setparam(D.groupelem{p},'stomesh',stomesh);
end

Ds = D;
for p=1:Ds.nbgroupelem
    
    elemin = getelem(Ds.groupelem{p},einglob,'global') ;
    elemcut = getelem(Ds.groupelem{p},ecutglob,'global') ;
    elemout = getelem(Ds.groupelem{p},eoutglob,'global') ;
    elemin = setlsenrich(elemin,0);
    elemout = setlsenrich(elemout,0);
    elemcut = setlsenrich(elemcut,1);
    elemin = setlstype(elemin,'in');
    
    elemin = setlsenrich(elemin,0);
    elemin = setlstype(elemin,'in');
    elemin = setlsnumber(elemin,getnumber(ls));
    enrich = ~ischarin('noenrich',varargin) & ~isempty(getmaterial(ls))  ;
    
    elemcut = setlsenrich(elemcut,enrich);
    elemcut = setlstype(elemcut,'cut');
    elemcut = setlsnumber(elemcut,getnumber(ls));
    if isempty(getmaterial(ls))
        elemout = setlstype(elemout,'out');
    end
    
    if isempty(getlsnumber(elemout))
        elemout = setlsnumber(elemout,getnumber(ls));
    end
    
    Ds.groupelem([p,end+1,end+2]) = {elemin,elemcut,elemout};
end
Ds.nbgroupelem = length(Ds.groupelem);

Ds=removeemptygroup(Ds);

function Iin = find_Iin(xi,stodim,U,X,lspcval,varargin)

x=[];
for i=1:stodim
    x = [x transfer(U{i},X{i},xi(:,i))];
end

Iin=sign(full(double(randomeval(lspcval,x,X))));
%keyboard
%    Iin=full(double(randomeval(lspcval,x,X)));
%    keyboard
%    h = max(Iin,[],1)-min(Iin,[],1);
%    h = max(h);
%    rep=find(abs(Iin(:))/h<1e-2);
%    Iin(rep)=0;
%    Iin = sign(Iin);
return

function I_elemtot = recursivecut_ND(I_elemtot,I_elem,tol,varargin)
stodim=getcharin('stodim',varargin);

rep = comp_I_ND(I_elem.signsommet,stodim);

if size(I_elem.order,1)==0
    crit = 0;
else
    crit = (max(I_elem.order(end,:)));
end
if ((1/2)^(crit)<=tol) || rep==1
    I = I_elem.signsommet;
    if all(sum((I==1 | I==0),1) == size(I,1))
        I_elem.state = 1;
    elseif all(sum((I==-1 | I==0),1) == size(I,1))
        I_elem.state = -1;
    else
        I_elem.state = 0;
    end
    
    I_elemtot = [I_elemtot , {I_elem}];
else
    I_subelem = cut_I_elem_ani_ND(I_elem,ones(1,stodim),1:stodim,varargin{:});
    for i=1:length(I_subelem)
        I_elemtot=recursivecut_ND(I_elemtot,I_subelem{i},tol,varargin{:});
    end
end

return

function I_elemtot = recursivecut_ani_ND(I_elemtot,I_elem,tol,varargin)
stodim=getcharin('stodim',varargin);
if ischarin('ani+',varargin)
    ani=2;
else
    ani=1;
end
if ani==1
    [cut cut1]= comp_I_ani_ND(I_elem.signsommet,stodim,tol,I_elem.order);
elseif ani==2
    [cut cut1]= comp_I_ani_ND_2(I_elem.signsommet,stodim,tol,I_elem.order);
end

if size(I_elem.order,1)==0
    crit = 0;
else
    crit = I_elem.order(end,:);
end
if all(((1/2).^crit)<=tol) || isempty(cut1)
    I = I_elem.signsommet;
    if all(sum((I==-1 | I==0),1) == size(I,1))
        I_elem.state = -1;
    elseif all(sum((I==1 | I==0),1) == size(I,1))
        I_elem.state = 1;
    else
        I_elem.state = 0;
    end
    
    I_elemtot = [I_elemtot , {I_elem}];
else
    I_subelem = cut_I_elem_ani_ND(I_elem,cut,cut1,varargin{:});
    for i=1:length(I_subelem)
        I_elemtot=recursivecut_ani_ND(I_elemtot,I_subelem{i},tol,varargin{:});
    end
end

return

function I_subelem = cut_I_elem_ani_ND(I_elem,cut,cut1,varargin)

stodim=getcharin('stodim',varargin);
varargin = delcharin('stodim',varargin);
nbsubelem = 2^(size(cut1,2));
I = calc_multi_indices(stodim,1,2);
I=I(:,1:end-1);
way = calc_way(I,cut,nbsubelem,stodim);
for i=1:nbsubelem
    if size(I_elem.order,1)==0
        I_subelem{i}.order = cut;
    else
        I_subelem{i}.order = [I_elem.order;I_elem.order(end,:)+cut];
    end
    I_subelem{i}.way = [I_elem.way;way(i,:)];
end
xisub = zeros(0,stodim);
for j=1:length(I_subelem)
    xisub=[xisub;calc_sommets(calc_xi_ani(I_subelem{j},stodim),stodim)];
end
Iin = find_Iin(xisub,stodim,varargin{:});

for j=1:length(I_subelem)
    I_subelem{j}.signsommet = Iin(:,(j-1)*2^stodim+[1:2^stodim]);
end


return

function rep = comp_I_ND(Iin,stodim)

rep = 1;
for i=2:2^stodim
    rep = rep & all(Iin(:,i)==Iin(:,1));
end

return

function [rep rep1] = comp_I_ani_ND(Iin,stodim,tol,order)
% 0 on ne decoupe pas 1 on decoupe
if size(order,1)==0
    order = zeros(1,size(order,2));
else
end

rep = zeros(1,stodim);
rep1 = [];
I=calc_multi_indices(stodim,1,2);
I=I(:,1:end-1);
for i=1:stodim
    [Ii,repi]=sortrows(I,i);
    faces{i}=reshape(repi,2^(stodim-1),2);
end

for i=1:stodim
    F1 = faces{i}(:,1);
    F2 = faces{i}(:,2);
    if any(any((Iin(:,F1)==Iin(:,F2))==0)) && (1/2).^(order(end,i))>tol
        rep(1,i)=1;
        rep1 = [rep1 i];
    end
end
return

function [rep rep1] = comp_I_ani_ND_2(Iin,stodim,tol,order)
% 0 on ne decoupe pas 1 on decoupe
if size(order,1)==0
    order = zeros(1,size(order,2));
else
end
rep = zeros(1,stodim);
rep1 = [];
if all(Iin(:,1)==Iin(:,2))
    [resp1] = calc_Iin(Iin(:,1),Iin,stodim);
    if resp1==0
    elseif resp1 ==1
        decoup = test_tol(order,stodim,tol);
        rep(1,decoup)=1;
        rep1 = [rep1 decoup];
    elseif resp1==2
        I=calc_multi_indices(stodim,1,2);
        I=I(:,1:end-1);
        for i=1:stodim
            [Ii,repi]=sortrows(I,i);
            faces{i}=reshape(repi,2^(stodim-1),2);
        end
        
        for i=1:stodim
            F1 = faces{i}(:,1);
            F2 = faces{i}(:,2);
            if any(any((Iin(:,F1)==Iin(:,F2))==0)) && (1/2).^(order(end,i))>tol
                rep(1,i)=1;
                rep1 = [rep1 i];
            end
        end
    end
else
    [resp1] = calc_Iin(Iin(:,1),Iin,stodim);
    [resp2] = calc_Iin(Iin(:,2),Iin,stodim);
    if resp1==0 || resp2==0
        decoup = test_tol(order,stodim,tol);
        rep(1,decoup)=1;
        rep1 = [rep1 decoup];
    else
        I=calc_multi_indices(stodim,1,2);
        I=I(:,1:end-1);
        for i=1:stodim
            [Ii,repi]=sortrows(I,i);
            faces{i}=reshape(repi,2^(stodim-1),2);
        end
        
        for i=1:stodim
            F1 = faces{i}(:,1);
            F2 = faces{i}(:,2);
            if any(any((Iin(:,F1)==Iin(:,F2))==0)) && (1/2).^(order(end,i))>tol
                rep(1,i)=1;
                rep1 = [rep1 i];
            end
        end
    end
end
return

function [xi] = calc_xi_ani(I_elem,stodim)

xi = sum((I_elem.way-1).*((1/2).^(I_elem.order)),1);
if size(I_elem.order,1)==0
    xi = [xi;xi+(1/2).^zeros(1,stodim)];
else
    xi = [xi;xi+(1/2).^(I_elem.order(end,:))];
end

return

function [xi] = calc_sommets(xiuni,stodim)
I=calc_multi_indices(stodim,1,2);
I=I(:,1:end-1);
xi=[];
for i=1:stodim
    xi = [xi,xiuni(I(:,i)+1,i)];
end

return

function rep = I_elemtotin(I_elemtot)
state=sparse(1,length(I_elemtot));
for i=1:length(I_elemtot)
    state(i)= I_elemtot{i}.state;
end
if all(state==-1)
    rep =1;
else
    rep=0;
end
return

function rep =I_elemtotout(I_elemtot)
state=sparse(1,length(I_elemtot));
for i=1:length(I_elemtot)
    state(i)= I_elemtot{i}.state;
end
if all(state==1)
    rep =1;
else
    rep=0;
end
return

function way = calc_way(I,cut,nbsubelem,stodim)
way = sortrows(I,find(cut==0));
way = way(1:nbsubelem,:);
for i=1:stodim
    if cut(i)==1
        way(:,i) = way(:,i)+1;
    elseif cut(i)==0
        way(:,i) = ones(size(way,1),1);
    else
        error('c''est 1 ou 0 !')
    end
end

return

function [rep1 rep] = calc_Iin(Iintest,Iin,stodim)

rep=[];
for i =3:2^(stodim)
    rep = [rep ;all(Iintest==Iin(:,i))];
end
if sum(rep)==2^(stodim)-2
    rep1=0; % on ne decoupe pas
elseif sum(rep)==2^(stodim)-3
    rep1=1; % on decoupe dans une seule direction
else
    rep1=2; % on decoupe mormalement
end
return

function rep = test_tol(order,stodim,tol)
rep = zeros(1,stodim);
for j=1:stodim
    rep(j) = ((1/2).^(order(end,j))>tol)*j;
end
rep = min(find(rep));
return
