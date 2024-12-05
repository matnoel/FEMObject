%% Homogeneisation meca 2D separee en XY, materiaux aleatoires
% Definition des modeles
dim=2;
size_img=[1 1];
nelem=[30 31];

Dx=DOMAIN(1,0,size_img(1));
Dy=DOMAIN(1,0,size_img(2));
V=getvolume(Dx)*getvolume(Dy);

Sx=mesh(Dx,nelem(1));
Sy=mesh(Dy,nelem(2));

Sx=createddlnode(Sx,DDL('u'));
Sy=createddlnode(Sy,DDL('u'));

%% Creation de l'indicatrice
[X,Y]=meshgrid( linspace(0,size_img(1),nelem(1)+1),...
    linspace(0,size_img(2),nelem(2)+1));

cx=0.5;
cy=0.5;
rfibre=0.3338;

lste=sqrt((X-cx).^2+(Y-cy).^2)-rfibre;

delta=2/50;
indicf=1-(1+tanh(lste/delta))/2;

indicfs=multisvd(indicf,'maxorder',50,'tol',1e-3);
indicfs=permutedim(indicfs,[2;1]); % X et Y sont remis dans l'ordre

%% Introduction des variables aleatoires
C1=[1 0 0;0 1 0; 0 0 0.5];
C2=[0 1 0;1 0 0; 0 0 -0.5];

Ef=10;
rvEf=RVNORMAL(Ef,Ef/10);
nuf=0.3;
nf=((1+nuf)*(1-2*nuf));
cf11=(1-nuf)/nf;
cf12=nuf/nf;
cf=cf11*C1+cf12*C2;

Em=1;
rvEm=RVNORMAL(Em,Em/10);
num=0.3;
nm=((1+num)*(1-2*num));
cm11=(1-num)/nm;
cm12=num/nm;
cm=cm11*C1+cm12*C2;

R=RANDVARS(rvEf,rvEm);

[X,PC]=PCTPMODEL(R,'order',6,'pcg');
% [X,PC]=PCTPMODEL(R,'order',6,'pcg','groups',{1:2;3:4});

oPC=one(PC);
oS=SEPMATRIX(oPC);
o=getphi(oPC);
for i=1:numel(o)
    o{i}=full(o{i}); %Pbs avec le format sparse
end

ef = calc_ximasse(X{1});
ef = get_ximasse(ef);
pf = getphi(X{1});
em=calc_ximasse(X{2});
em=get_ximasse(em);
pm=getphi(X{2});

%% Construction de l'operateur
Ex=zeros(3,2);
Ex(1,1)=1;
Ex(3,2)=1;

Ey=zeros(3,2);
Ey(2,2)=1;
Ey(3,1)=1;

kx=setfree(BILINFORM(1,1,1),0);
kx=kx{Sx}(:,:);
mx=setfree(BILINFORM(0,0,1),0);
mx=mx{Sx}(:,:);
dxl=setfree(BILINFORM(1,0,1),0);
dxl=dxl{Sx}(:,:);
dxr=dxl';

ky=setfree(BILINFORM(1,1,1),0);
ky=ky{Sy}(:,:);
my=setfree(BILINFORM(0,0,1),0);
my=my{Sy}(:,:);
dyl=setfree(BILINFORM(1,0,1),0);
dyl=dyl{Sy}(:,:);
dyr=dyl';

Gxxm=Ex'*cm*Ex;
Gyym=Ey'*cm*Ey;
Gxym=Ex'*cm*Ey;
Gyxm=Gxym';

Gxxf=Ex'*cf*Ex;
Gyyf=Ey'*cf*Ey;
Gxyf=Ex'*cf*Ey;
Gyxf=Gxyf';

A1=SEPMATRIX({kx,my,Gxxm,em{:};...
    mx,ky,Gyym,em{:};...
    dxl,dyr,Gxym,em{:};...
    dxr,dyl,Gyxm,em{:}});

%As,A1s,A2s,A3s sont utilises plus tard pour
% faire des tirages aleatoires de l'operateur

A1s=SEPMATRIX({kx,my,Gxxm,pm{:};...
    mx,ky,Gyym,pm{:};...
    dxl,dyr,Gxym,pm{:};...
    dxr,dyl,Gyxm,pm{:}});


A2=SEPMATRIX(getdim(A1));
A2s=SEPMATRIX(getdim(A1));
A3=SEPMATRIX(getdim(A1));
A3s=SEPMATRIX(getdim(A1));
for i=1:indicfs.m
    kx=setfree(BILINFORM(1,1,indicfs.F{i,1},0),0);
    kx=kx{Sx}(:,:);
    mx=setfree(BILINFORM(0,0,indicfs.F{i,1},0),0);
    mx=mx{Sx}(:,:);
    dxl=setfree(BILINFORM(1,0,indicfs.F{i,1},0),0);
    dxl=dxl{Sx}(:,:);
    dxr=dxl';

    ky=setfree(BILINFORM(1,1,indicfs.F{i,2},0),0);
    ky=ky{Sy}(:,:);
    my=setfree(BILINFORM(0,0,indicfs.F{i,2},0),0);
    my=my{Sy}(:,:);
    dyl=setfree(BILINFORM(1,0,indicfs.F{i,2},0),0);
    dyl=dyl{Sy}(:,:);
    dyr=dyl';

    A2=A2+SEPMATRIX({ kx,my,Gxxf,ef{:};...
        mx,ky,Gyyf,ef{:};...
        dxl,dyr,Gxyf,ef{:};...
        dxr,dyl,Gyxf,ef{:}},...
        [indicfs.alpha(i),indicfs.alpha(i),indicfs.alpha(i),indicfs.alpha(i)]);
    A2s=A2s+SEPMATRIX({ kx,my,Gxxf,pf{:};...
        mx,ky,Gyyf,pf{:};...
        dxl,dyr,Gxyf,pf{:};...
        dxr,dyl,Gyxf,pf{:}},...
        [indicfs.alpha(i),indicfs.alpha(i),indicfs.alpha(i),indicfs.alpha(i)]);

    A3=A3+SEPMATRIX({ kx,my,Gxxm,em{:};...
        mx,ky,Gyym,em{:};...
        dxl,dyr,Gxym,em{:};...
        dxr,dyl,Gyxm,em{:}},...
        [indicfs.alpha(i),indicfs.alpha(i),indicfs.alpha(i),indicfs.alpha(i)]);
    A3s=A3s+SEPMATRIX({ kx,my,Gxxm,pm{:};...
        mx,ky,Gyym,pm{:};...
        dxl,dyr,Gxym,pm{:};...
        dxr,dyl,Gyxm,pm{:}},...
        [indicfs.alpha(i),indicfs.alpha(i),indicfs.alpha(i),indicfs.alpha(i)]);
end
A=A1+A2-A3;
As=A1s+A2s-A3s;
As=PCTPMATRIX(As,PC,4:5,'noexpand',1);

%% Conditions aux limites de Dirichlet
Sxe=Sx;
Sye=Sy;

Sxe=addcl(Sxe,POINT([0;size_img(1)]),'u');
Sye=addcl(Sye,POINT([0;size_img(2)]),'u');

Af=freematrix(A,1,Sxe);
Af=freematrix(Af,2,Sye);

x = getcoord(getnode(Sxe));
y = getcoord(getnode(Sye));

% Mode 1
ubcx = SEPMATRIX({x(:),ones(length(y),1),[1;0],o{:}});
b=-A*ubcx;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,...
    'update',1,'updatedim',3:getdim(A));
[u1sep,result1ref] = solve(Af,b,PGD);
u1e=unfreevector(u1sep,1,Sxe);
u1e=unfreevector(u1e,2,Sye);
u1e=ubcx+u1e;
u1e=PCTPMATRIX(u1e,PC,4:5,'noexpand',1);

% Mode 2
ubcy = SEPMATRIX({ones(length(x),1),y(:),[0;1],o{:}});
b=-A*ubcy;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,...
    'update',1,'updatedim',3:getdim(A));
[u2sep,result2ref] = solve(Af,b,PGD);
u2e=unfreevector(u2sep,1,Sxe);
u2e=unfreevector(u2e,2,Sye);
u2e=ubcy+u2e;
u2e=PCTPMATRIX(u2e,PC,4:5,'noexpand',1);

% Mode 3
ubcm = SEPMATRIX({x(:),ones(length(y),1),[0;1],o{:}});
b=-A*ubcm;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,...
    'update',1,'updatedim',3:getdim(A));
[u3sep,result3ref] = solve(Af,b,PGD);
u3e=unfreevector(u3sep,1,Sxe);
u3e=unfreevector(u3e,2,Sye);
u3e=ubcm+u3e;
u3e=PCTPMATRIX(u3e,PC,4:5,'noexpand',1);

%% Conditions aux limites periodiques
Sxp=Sx;
Syp=Sy;

Sxp=addclperiodic(Sxp,POINT(0),POINT(size_img(1)),'u');
Syp=addclperiodic(Syp,POINT(0),POINT(size_img(2)),'u');

% Blocage des translations de solide rigide
lx = setfree( LINFORM(0,1),0);
Lx = calc_vector(lx,Sxp);
Ly = calc_vector(lx,Syp);
alpha=1;
ff=oS.F(2:oS.dim);
Ap = A + SEPMATRIX({Lx*Lx',Ly*Ly',eye(2),ff{:}},alpha);

Af=freematrix(Ap,1,Sxp);
Af=freematrix(Af,2,Syp);

% Mode 1
b=-A*ubcx;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,...
    'update',1,'updatedim',3:getdim(A));
[u4sep,result4ref] = solve(Af,b,PGD);
u1p=unfreevector(u4sep,1,Sxp);
u1p=unfreevector(u1p,2,Syp);
u1p=ubcx+u1p;
u1p=PCTPMATRIX(u1p,PC,4:5,'noexpand',1);

% Mode 2
b=-A*ubcy;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,...
    'update',1,'updatedim',3:getdim(A));
[u5sep,result5ref] = solve(Af,b,PGD);
u2p=unfreevector(u5sep,1,Sxp);
u2p=unfreevector(u2p,2,Syp);
u2p=ubcy+u2p;
u2p=PCTPMATRIX(u2p,PC,4:5,'noexpand',1);

% Mode 3
b=-A*ubcm;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,...
    'update',1,'updatedim',3:getdim(A));
[u6sep,result6ref] = solve(Af,b,PGD);
u3p=unfreevector(u6sep,1,Sxp);
u3p=unfreevector(u3p,2,Syp);
u3p=ubcm+u3p;
u3p=PCTPMATRIX(u3p,PC,4:5,'noexpand',1);

%% Conditions aux limites de Neumann
Sxn=Sx;
Syn=Sy;

% Blocage des rotations de solide rigide
x=getcoord(getnode(Sxn));
y=getcoord(getnode(Syn));

mx = setfree( BILINFORM(0,0,1),0);
Mx = calc_matrix(mx,Sxn);
My = calc_matrix(mx,Syn);

d11 = [1 0;0 0];
d12 = [0 1;0 0];
d21 = [0 0;1 0];
d22 = [0 0;0 1];

App=Ap+ ...
    SEPMATRIX({Lx*Lx', My*(y*y')*My',d11,ff{:};...
Lx*x'*Mx',My*y*Ly',d12,ff{:};...
Mx*x*Lx',Ly*y'*My',d21,ff{:};...
Mx*(x*x')*Mx',Ly*Ly',d22,ff{:}},[1 -1 -1 1]);

l = setfree(LINFORM(0,1),0);
% Mode 1
ly=calc_vector(l,Syn);
xm=zeros(nelem(1)+1,1);
xm(1)=-1;
xp=zeros(nelem(1)+1,1);
xp(nelem(1)+1)=1;
b = SEPMATRIX({xm,ly,[1;0],o{:};xp,ly,[1;0],o{:}});
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,...
    'update',1,'updatedim',3:getdim(A));
[u7sep,result7ref] = solve(App,b,PGD);
u1n=u7sep;
u1n=PCTPMATRIX(u1n,PC,4:5,'noexpand',1);

% Mode 2
lx=calc_vector(l,Sxn);
ym=zeros(nelem(2)+1,1);
ym(1)=-1;
yp=zeros(nelem(2)+1,1);
yp(nelem(2)+1)=1;
b = SEPMATRIX({lx,ym,[0;1],o{:};lx,yp,[0;1],o{:}});
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,...
    'update',1,'updatedim',3:getdim(A));
[u8sep,result8ref] = solve(App,b,PGD);
u2n=u8sep;
u2n=PCTPMATRIX(u2n,PC,4:5,'noexpand',1);

% Mode 3
b = SEPMATRIX({xm,ly,[0;1],o{:};xp,ly,[0;1],o{:};...
    lx,ym,[1;0],o{:};lx,yp,[1;0],o{:}});
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,...
    'update',1,'updatedim',3:getdim(A));
[u9sep,result9ref] = solve(App,b,PGD);
u3n=u9sep;
u3n=PCTPMATRIX(u3n,PC,4:5,'noexpand',1);

%% Determination du tenseur d'elasticite homogeneise
% Utilisation de l'energie de deformation

nbreal=10;
x=zeros(nbreal,2);
Emat=zeros(nbreal,1);
Efib=zeros(nbreal,1);

Cer=zeros(nbreal,9);
Cpr=zeros(nbreal,9);
Cnr=zeros(nbreal,9);

cc=zeros(3);

fprintf('Calcul des conductivites\n')

for i=1:nbreal
    [Efib(i),x(i,:)]=random(X{1});
    Emat(i)=randomeval(X{2},x(i,:));

    aa=randomeval(As,x(i,:));

    u={randomeval(u1e,x(i,:)),...
        randomeval(u2e,x(i,:)),... 
        randomeval(u3e,x(i,:))};
    v={expand(aa*u{1}),...
        expand(aa*u{2}),...
        expand(aa*u{3})};
    for j=1:3
        u{j}=expand(u{j});
    end
    for j=1:3
        for k=1:3
            cc(j,k)=u{j}(:)'*v{k}(:);
        end
    end
    Cer(i,:)=cc(:);

    u={randomeval(u1p,x(i,:)),...
        randomeval(u2p,x(i,:)),... 
        randomeval(u3p,x(i,:))};
    v={expand(aa*u{1}),...
        expand(aa*u{2}),...
        expand(aa*u{3})};
    for j=1:3
        u{j}=expand(u{j});
    end
    for j=1:3
        for k=1:3
            cc(j,k)=u{j}(:)'*v{k}(:);
        end
    end
    Cpr(i,:)=cc(:);

    u={randomeval(u1n,x(i,:)),...
        randomeval(u2n,x(i,:)),... 
        randomeval(u3n,x(i,:))};
    v={expand(aa*u{1}),...
        expand(aa*u{2}),...
        expand(aa*u{3})};
    for j=1:3
        u{j}=expand(u{j});
    end
    for j=1:3
        for k=1:3
            cc(j,k)=u{j}(:)'*v{k}(:);
        end
    end
    cc=inv(cc);
    Cnr(i,:)=cc(:);
    fprintf('%6.0f%%\n',100*i/nbreal)
end

Cer=Cer/V;
Cpr=Cpr/V;
Cnr=Cnr/V;

C=zeros(nbreal,3);
C(:,1)=Cer(:,1);
C(:,2)=Cpr(:,1);
C(:,3)=Cnr(:,1);

disp([Emat Efib C])

%% Controle de la premiere valeur
% Comparaison entre une solution deterministe et un tirage aleatoire
% effectue precedemment

Cf=elasticity_tensor(Efib(1),nuf,2);
Cm=elasticity_tensor(Emat(1),num,2);
Cd=Cf-Cm;

C1=SEPMATRIX(3);
C1.alpha=1;
C1.F{1,1}=ones(nelem(1)+1,1);
C1.F{1,2}=ones(nelem(2)+1,1);
C1.F{1,3}=Cm;

C2=SEPMATRIX(3);
C2.alpha=indicfs.alpha;
C2.F=indicfs.F;
for i=1:indicfs.m
    C2.F{i,3}=Cd;
end

Cmap=C1+C2;

Ex=zeros(3,2);
Ex(1,1)=1;
Ex(3,2)=1;

Ey=zeros(3,2);
Ey(2,2)=1;
Ey(3,1)=1;

Ad=SEPMATRIX(3); % Operateur deterministe
for i=1:Cmap.m
    kx=setfree(BILINFORM(1,1,Cmap.F{i,1},0),0);
    kx=kx{Sx}(:,:);
    mx=setfree(BILINFORM(0,0,Cmap.F{i,1},0),0);
    mx=mx{Sx}(:,:);
    dxl=setfree(BILINFORM(1,0,Cmap.F{i,1},0),0);
    dxl=dxl{Sx}(:,:);
    dxr=dxl';

    ky=setfree(BILINFORM(1,1,Cmap.F{i,2},0),0);
    ky=ky{Sy}(:,:);
    my=setfree(BILINFORM(0,0,Cmap.F{i,2},0),0);
    my=my{Sy}(:,:);
    dyl=setfree(BILINFORM(1,0,Cmap.F{i,2},0),0);
    dyl=dyl{Sy}(:,:);
    dyr=dyl';

    Gxx=Ex'*Cmap.F{i,3}*Ex;
    Gyy=Ey'*Cmap.F{i,3}*Ey;
    Gxy=Ex'*Cmap.F{i,3}*Ey;
    Gyx=Gxy';

    Ad=Ad+SEPMATRIX({kx,my,Gxx;mx,ky,Gyy;dxl,dyr,Gxy;dxr,dyl,Gyx},...
        [Cmap.alpha(i),Cmap.alpha(i),Cmap.alpha(i),Cmap.alpha(i)]);
end

Asd=randomeval(As,x(1,:)); % Tirage aleatoire de l'operateur

%% Conditions aux limites de Dirichlet
Sxe=Sx;
Sye=Sy;

Sxe=addcl(Sxe,POINT([0;size_img(1)]),'u');
Sye=addcl(Sye,POINT([0;size_img(2)]),'u');

Af=freematrix(Ad,1,Sxe);
Af=freematrix(Af,2,Sye);
Afs=freematrix(Asd,1,Sxe);
Afs=freematrix(Afs,2,Sye);

x = getcoord(getnode(Sxe));
y = getcoord(getnode(Sye));

% Mode 1
ubcx = SEPMATRIX({x(:),ones(length(y),1),[1;0]});
b=-Ad*ubcx;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(Ad),'maxorder',30,'tol',1e-5);
[u1sepd,result1refd] = solve(Af,b,PGD);
u1ed=unfreevector(u1sepd,1,Sxe);
u1ed=unfreevector(u1ed,2,Sye);
u1ed=ubcx+u1ed; % Solution deterministe

[u1sepsd,result1refsd] = solve(Afs,b,PGD);
u1esd=unfreevector(u1sepsd,1,Sxe);
u1esd=unfreevector(u1esd,2,Sye);
u1esd=ubcx+u1esd; % Solution provenant du tirage aleatoire

ud=expand(u1ed);
usd=expand(u1esd);
e1e=norm(ud(:)-usd(:),'fro')/norm(ud(:),'fro'); % Erreur entre les solutions

% Mode 2
ubcy = SEPMATRIX({ones(length(x),1),y(:),[0;1]});
b=-Ad*ubcy;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(Ad),'maxorder',30,'tol',1e-5);
[u2sepd,result2refd] = solve(Af,b,PGD);
u2ed=unfreevector(u2sepd,1,Sxe);
u2ed=unfreevector(u2ed,2,Sye);
u2ed=ubcy+u2ed;

[u2sepsd,result2refsd] = solve(Afs,b,PGD);
u2esd=unfreevector(u2sepsd,1,Sxe);
u2esd=unfreevector(u2esd,2,Sye);
u2esd=ubcy+u2esd;

ud=expand(u2ed);
usd=expand(u2esd);
e2e=norm(ud(:)-usd(:),'fro')/norm(ud(:),'fro');

% Mode 3
ubcm = SEPMATRIX({x(:),ones(length(y),1),[0;1]});
b=-Ad*ubcm;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(Ad),'maxorder',30,'tol',1e-5);
[u3sepd,result3refd] = solve(Af,b,PGD);
u3ed=unfreevector(u3sepd,1,Sxe);
u3ed=unfreevector(u3ed,2,Sye);
u3ed=ubcm+u3ed;

[u3sepsd,result3refsd] = solve(Afs,b,PGD);
u3esd=unfreevector(u3sepsd,1,Sxe);
u3esd=unfreevector(u3esd,2,Sye);
u3esd=ubcm+u3esd;

ud=expand(u3ed);
usd=expand(u3esd);
e3e=norm(ud(:)-usd(:),'fro')/norm(ud(:),'fro');

%% Conditions aux limites periodiques
Sxp=Sx;
Syp=Sy;

Sxp=addclperiodic(Sxp,POINT(0),POINT(size_img(1)),'u');
Syp=addclperiodic(Syp,POINT(0),POINT(size_img(2)),'u');

lx = setfree( LINFORM(0,1),0);
Lx = calc_vector(lx,Sxp);
Ly = calc_vector(lx,Syp);
alpha=1;
Apd = Ad + SEPMATRIX({Lx*Lx',Ly*Ly',eye(2)},alpha);

Apsd = Asd + SEPMATRIX({Lx*Lx',Ly*Ly',eye(2)},alpha);

Af=freematrix(Apd,1,Sxp);
Af=freematrix(Af,2,Syp);

Afs=freematrix(Apsd,1,Sxp);
Afs=freematrix(Afs,2,Syp);

% Mode 1
b=-Ad*ubcx;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(Ad),'maxorder',30,'tol',1e-5);
[u4sepd,result4refd] = solve(Af,b,PGD);
u1pd=unfreevector(u4sepd,1,Sxp);
u1pd=unfreevector(u1pd,2,Syp);
u1pd=ubcx+u1pd;

[u4sepsd,result4refsd] = solve(Afs,b,PGD);
u1psd=unfreevector(u4sepsd,1,Sxp);
u1psd=unfreevector(u1psd,2,Syp);
u1psd=ubcx+u1psd; % Solution provenant du tirage aleatoire

ud=expand(u1pd);
usd=expand(u1psd);
e1p=norm(ud(:)-usd(:),'fro')/norm(ud(:),'fro'); % Erreur entre les solutions


% Mode 2
b=-Ad*ubcy;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(Ad),'maxorder',30,'tol',1e-5);
[u5sepd,result5refd] = solve(Af,b,PGD);
u2pd=unfreevector(u5sepd,1,Sxp);
u2pd=unfreevector(u2pd,2,Syp);
u2pd=ubcy+u2pd;

[u5sepsd,result5refsd] = solve(Afs,b,PGD);
u2psd=unfreevector(u5sepsd,1,Sxp);
u2psd=unfreevector(u2psd,2,Syp);
u2psd=ubcy+u2psd;

ud=expand(u2pd);
usd=expand(u2psd);
e2p=norm(ud(:)-usd(:),'fro')/norm(ud(:),'fro');

% Mode 3
b=-Ad*ubcm;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(Ad),'maxorder',30,'tol',1e-5);
[u6sepd,result6refd] = solve(Af,b,PGD);
u3pd=unfreevector(u6sepd,1,Sxp);
u3pd=unfreevector(u3pd,2,Syp);
u3pd=ubcm+u3pd;

[u6sepsd,result6refsd] = solve(Afs,b,PGD);
u3psd=unfreevector(u6sepsd,1,Sxp);
u3psd=unfreevector(u3psd,2,Syp);
u3psd=ubcm+u3psd;

ud=expand(u3pd);
usd=expand(u3psd);
e3p=norm(ud(:)-usd(:),'fro')/norm(ud(:),'fro');

%% Conditions aux limites de Neumann
Sxn=Sx;
Syn=Sy;

% Blocage des rotations de solide rigide
x=getcoord(getnode(Sxn));
y=getcoord(getnode(Syn));

mx = setfree( BILINFORM(0,0,1),0);
Mx = calc_matrix(mx,Sxn);
My = calc_matrix(mx,Syn);

d11 = [1 0;0 0];
d12 = [0 1;0 0];
d21 = [0 0;1 0];
d22 = [0 0;0 1];

Appd=Apd+ ...
    SEPMATRIX({Lx*Lx', My*(y*y')*My',d11;...
Lx*x'*Mx',My*y*Ly',d12;...
Mx*x*Lx',Ly*y'*My',d21;...
Mx*(x*x')*Mx',Ly*Ly',d22},[1 -1 -1 1]);

Appsd = Apsd + ...
    SEPMATRIX({Lx*Lx', My*(y*y')*My',d11;...
Lx*x'*Mx',My*y*Ly',d12;...
Mx*x*Lx',Ly*y'*My',d21;...
Mx*(x*x')*Mx',Ly*Ly',d22},[1 -1 -1 1]);

l = setfree(LINFORM(0,1),0);
% Mode 1
ly=calc_vector(l,Syn);
xm=zeros(nelem(1)+1,1);
xm(1)=-1;
xp=zeros(nelem(1)+1,1);
xp(nelem(1)+1)=1;
b = SEPMATRIX({xm,ly,[1;0];xp,ly,[1;0]});
PGD = SEPSOLVER(getdim(Ad),'maxorder',50,'tol',1e-5);
[u7sepd,result7refd] = solve(Appd,b,PGD);
u1nd=u7sepd;

[u7sepsd,result7refsd] = solve(Appsd,b,PGD);
u1nsd=u7sepsd;

ud=expand(u1nd);
usd=expand(u1nsd);
e1n=norm(ud(:)-usd(:),'fro')/norm(ud(:),'fro');

% Mode 2
lx=calc_vector(l,Sxn);
ym=zeros(nelem(2)+1,1);
ym(1)=-1;
yp=zeros(nelem(2)+1,1);
yp(nelem(2)+1)=1;
b = SEPMATRIX({lx,ym,[0;1];lx,yp,[0;1]});
PGD = SEPSOLVER(getdim(A),'maxorder',50,'tol',1e-5);
[u8sepd,result8refd] = solve(Appd,b,PGD);
u2nd=u8sepd;

[u8sepsd,result8refsd] = solve(Appsd,b,PGD);
u2nsd=u8sepsd;

ud=expand(u2nd);
usd=expand(u2nsd);
e2n=norm(ud(:)-usd(:),'fro')/norm(ud(:),'fro');

% Mode 3
b = SEPMATRIX({xm,ly,[0;1];xp,ly,[0;1];...
    lx,ym,[1;0];lx,yp,[1;0]});
PGD = SEPSOLVER(getdim(Ad),'maxorder',50,'tol',1e-5);
[u9sepd,result9refd] = solve(Appd,b,PGD);
u3nd=u9sepd;

[u9sepsd,result9refsd] = solve(Appsd,b,PGD);
u3nsd=u9sepsd;

ud=expand(u3nd);
usd=expand(u3nsd);
e3n=norm(ud(:)-usd(:),'fro')/norm(ud(:),'fro');

%% Determination du tenseur d'elasticite homogeneisees
% Utilisation de l'energie de deformation

Ces=reshape(Cer(1,:),3,3);
Cps=reshape(Cpr(1,:),3,3);
Cns=reshape(Cnr(1,:),3,3);

u={u1ed,u2ed,u3ed};
Ced=zeros(3);
for i=1:3
    for j=1:3
        Ced(i,j)=expand(u{i}'*Ad*u{j});
    end
end
Ced=Ced/V;

u={u1esd,u2esd,u3esd};
Cesd=zeros(3);
for i=1:3
    for j=1:3
        Cesd(i,j)=expand(u{i}'*Asd*u{j});
    end
end
Cesd=Cesd/V;

u={u1pd,u2pd,u3pd};
Cpd=zeros(3);
for i=1:3
    for j=1:3
        Cpd(i,j)=expand(u{i}'*Ad*u{j});
    end
end
Cpd=Cpd/V;

u={u1psd,u2psd,u3psd};
Cpsd=zeros(3);
for i=1:3
    for j=1:3
        Cpsd(i,j)=expand(u{i}'*Asd*u{j});
    end
end
Cpsd=Cpsd/V;

u={u1nd,u2nd,u3nd};
Cnd=zeros(3);
for i=1:3
    for j=1:3
        Cnd(i,j)=expand(u{i}'*Ad*u{j});
    end
end
Cnd=Cnd/V;
Cnd=inv(Cnd);

u={u1nsd,u2nsd,u3nsd};
Cnsd=zeros(3);
for i=1:3
    for j=1:3
        Cnsd(i,j)=expand(u{i}'*Asd*u{j});
    end
end
Cnsd=Cnsd/V;
Cnsd=inv(Cnsd);

eCe=norm(Ced-Cesd,'fro')/norm(Ced,'fro');
eCp=norm(Cpd-Cpsd,'fro')/norm(Cpd,'fro');
eCn=norm(Cnd-Cnsd,'fro')/norm(Cnd,'fro');

eCes=norm(Ces-Cesd,'fro')/norm(Ces,'fro');
eCps=norm(Cps-Cpsd,'fro')/norm(Cps,'fro');
eCns=norm(Cns-Cnsd,'fro')/norm(Cns,'fro');


