function [Ds,randh] = lsrandomsplit(ls,tol,varargin)
% function [Ds,randh] = lsrandomsplit(ls,tol,varargin)

X = RANDVARS(ls);

if length(X)>1
    error('lssplitelem non programme en 2 dimensions stochastiques et plus')
end

D=getclassin('MODEL',varargin);
if isempty(D)
    D=ls.D;
end

ecutglob = [];
einglob = [];
eoutglob = [];

xitot=[0,1];

for p=1:D.nbgroupelem
    elem = D.groupelem{p};
    
    elem1 = repelemcutxi(0,ls,D,elem,X);
    elem2 = repelemcutxi(1,ls,D,elem,X);
    elemsto{1}=union(elem1{1},elem2{1});
    elemsto{2}=union(elem1{2},elem2{2});
    elemsto{3}=union(elem1{3},elem2{3});
    fprintf('Decomposition stochastique, tolerance %d ... ',tol)
    compteur = 1;
    if ~elemcmp(elem1,elem2)
        [xitot,elemsto,compteur]=recursivecut(compteur,0,1,tol,xitot,elemsto,elem1,elem2,ls,D,elem,X);
    end
    fprintf('\n')
    ecutglob = union(ecutglob, getnumber(elem,elemsto{2}));
    einglob = union(einglob, getnumber(elem,elemsto{1}));
    eoutglob = union(eoutglob, getnumber(elem,elemsto{3}));
    
end

xitot=sort(xitot);

separetol=max(10*tol,1e-10);
xitot;
repelim = find(xitot(3:end-1)-xitot(2:end-2)<separetol);
xitot(repelim+1) = [];
if length(xitot)>2 && (xitot(2)-xitot(1)<separetol)
    xitot(2)=[];
end
if length(xitot)>2 && (xitot(end)-xitot(end-1)<separetol)
    xitot(end-1)=[];
end
randh = POLYFE(xitot);

einglob=setdiff(einglob,ecutglob);
eoutglob=setdiff(eoutglob,ecutglob);

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

%figure(3)
%clf
%plot(Ds,'lsenrich','edges')

function [xitot,elemsto,compteur]=recursivecut(compteur,xi1,xi2,tol,xitot,elemsto,elem1,elem2,varargin)

xi=[xi1,(xi1+xi2)/2,xi2];

elem3 = repelemcutxi(xi(2),varargin{:});
elemsto{1} = union(elemsto{1},elem3{1});
elemsto{2} = union(elemsto{2},elem3{2});
elemsto{3} = union(elemsto{3},elem3{3});

if elemcmp(elem1,elem3) && elemcmp(elem2,elem3)
    xitot = [xitot,xi1,xi2];
else
    if (xi(2)-xi(1))<tol
        xitot = [xitot,xi(2)];
    else
        
        compteur = compteur+1;
        pourcentage(log(compteur+1)/log(2),ceil(log(tol)/log(1/2)));
        
        if ~elemcmp(elem1,elem3)
            [xitot,elemsto,compteur] = recursivecut(compteur,xi(1),xi(2),tol,xitot,elemsto,elem1,elem3,varargin{:});
        end
        if ~elemcmp(elem2,elem3)
            [xitot,elemsto,compteur] = recursivecut(compteur,xi(2),xi(3),tol,xitot,elemsto,elem3,elem2,varargin{:});
        end
    end
end

return

function e = repelemcutxi(xi,ls,D,elem,X)
% function e = repelemcutxi(xi,ls,D,elem,X)

XU = setnumber(RANDVARS(RVUNIFORM(0,1)),getnumber(X));

x=transfer(XU,X,xi);
ls = randomeval(ls,x,X);
ls=lseval(ls,D);
e{1} = find(lsisin(elem,ls,D.node));
e{2} = find(lsiscut(elem,ls,D.node));
e{3} = find(lsisout(elem,ls,D.node));


return

function ok = elemcmp(elem1,elem2)
% function ok = elemcmp(elem1,elem2)

if (length(elem1{1})==length(elem2{1})) && (length(elem1{2})==length(elem2{2})) && (length(elem1{3})==length(elem2{3}))
    ok = all(sort(elem1{1})==sort(elem2{1})) & all(sort(elem1{2})==sort(elem2{2})) & all(sort(elem1{3})==sort(elem2{3}));
else
    ok=0;
end

