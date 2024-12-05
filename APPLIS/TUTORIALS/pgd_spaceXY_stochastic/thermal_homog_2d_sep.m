%% Homogeneisation thermique 2D separee en XY
% Definition des modeles
dim=2;
size_img=[1 1];
nelem=[50 60];

Dx=DOMAIN(1,0,size_img(1));
Dy=DOMAIN(1,0,size_img(2));
V=getvolume(Dx)*getvolume(Dy);

Sx=mesh(Dx,nelem(1));
Sy=mesh(Dy,nelem(2));

Sx=createddlnode(Sx,DDL('T'));
Sy=createddlnode(Sy,DDL('T'));

%% Creation de l'indicatrice
D=DOMAIN(dim,[0 0],size_img);
V=getvolume(D);
S=mesh(D,nelem(1),nelem(2));
S=createddlnode(S,DDL('T'));

cx=0.5;
cy=0.5;
rfibre=0.3338;

lst=LSCIRCLE(cx,cy,rfibre);
lste=double(lseval(lst,S));

delta=2/50;
indicf1=1-(1+tanh(lste/delta))/2;
indicf1=reshape(indicf1,nelem(2)+1,nelem(1)+1);
%indicf=double(lste<0);

%% Creation de l'indicatrice V2 
% (bcp plus rapide pour les gros modeles)
% ATTENTION: inversion des dimensions (X->Y)
[X,Y]=meshgrid( linspace(0,size_img(1),nelem(1)+1),...
    linspace(0,size_img(2),nelem(2)+1));

lste=sqrt((X-cx).^2+(Y-cy).^2)-rfibre;

delta=2/50;
indicf2=1-(1+tanh(lste/delta))/2;

%% Separation de l'indicatrice et creation de la carte de materiaux
indicfs=multisvd(indicf2,'maxorder',50,'tol',1e-3);

sone=SEPMATRIX(2);
sone.alpha=1;
sone.F{1,1}=ones(nelem(2)+1,1);
sone.F{1,2}=ones(nelem(1)+1,1);

kf=10;
km=1;
kmap=km*sone+(kf-km)*indicfs;

kmap=permutedim(kmap,[2;1]); % X et Y sont remis dans l'ordre

%% Construction de l'operateur
A=SEPMATRIX(2);
for i=1:kmap.m
    kx=setfree(BILINFORM(1,1,kmap.F{i,1},0),0);
    kx=kx{Sx}(:,:);
    ky=setfree(BILINFORM(1,1,kmap.F{i,2},0),0);
    ky=ky{Sy}(:,:);
    mx=setfree(BILINFORM(0,0,kmap.F{i,1},0),0);
    mx=mx{Sx}(:,:);
    my=setfree(BILINFORM(0,0,kmap.F{i,2},0),0);
    my=my{Sy}(:,:);
    A=A+SEPMATRIX({kx,my;mx,ky},[kmap.alpha(i),kmap.alpha(i)]);
end

%% Conditions aux limites de Dirichlet
Sxe=Sx;
Sye=Sy;

Sxe=addcl(Sxe,POINT([0;size_img(1)]),'T');
Sye=addcl(Sye,POINT([0;size_img(2)]),'T');

Af=freematrix(A,1,Sxe);
Af=freematrix(Af,2,Sye);

x = getcoord(getnode(Sxe));
y = getcoord(getnode(Sye));

% Mode 1
ubcx = SEPMATRIX({x(:),ones(length(y),1)});
b=-A*ubcx;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u1sep,result1ref] = solve(Af,b,PGD);
u1e=unfreevector(u1sep,1,Sxe);
u1e=unfreevector(u1e,2,Sye);
u1e=ubcx+u1e;

% Mode 2
ubcy = SEPMATRIX({ones(length(x),1),y(:)});
b=-A*ubcy;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u2sep,result2ref] = solve(Af,b,PGD);
u2e=unfreevector(u2sep,1,Sxe);
u2e=unfreevector(u2e,2,Sye);
u2e=ubcy+u2e;

%% Conditions aux limites periodiques
Sxp=Sx;
Syp=Sy;

Sxp=addclperiodic(Sxp,POINT(0),POINT(size_img(1)),'T');
Syp=addclperiodic(Syp,POINT(0),POINT(size_img(2)),'T');

lx = setfree( LINFORM(0,1),0);
Lx = calc_vector(lx,Sxp);
Ly = calc_vector(lx,Syp);
alpha=1;
Ap = A + SEPMATRIX({Lx*Lx',Ly*Ly'},alpha);

Af=freematrix(Ap,1,Sxp);
Af=freematrix(Af,2,Syp);

% Mode 1
b=-A*ubcx;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u3sep,result3ref] = solve(Af,b,PGD);
u1p=unfreevector(u3sep,1,Sxp);
u1p=unfreevector(u1p,2,Syp);
u1p=ubcx+u1p;

% Mode 2
b=-A*ubcy;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u4sep,result4ref] = solve(Af,b,PGD);
u2p=unfreevector(u4sep,1,Sxp);
u2p=unfreevector(u2p,2,Syp);
u2p=ubcy+u2p;

%% Conditions aux limites de Neumann
Sxn=Sx;
Syn=Sy;

l = setfree(LINFORM(0,1),0);
% Mode 1
ly=calc_vector(l,Syn);
xm=zeros(nelem(1)+1,1);
xm(1)=-1;
xp=zeros(nelem(1)+1,1);
xp(nelem(1)+1)=1;
b = SEPMATRIX({xm,ly;xp,ly});
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u5sep,result5ref] = solve(Ap,b,PGD);
u1n=u5sep;

% Mode 2
lx=calc_vector(l,Sxn);
ym=zeros(nelem(2)+1,1);
ym(1)=-1;
yp=zeros(nelem(2)+1,1);
yp(nelem(2)+1)=1;
b = SEPMATRIX({lx,ym;lx,yp});
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u6sep,result6ref] = solve(Ap,b,PGD);
u2n=u6sep;

%% Determination des conductivites homogeneisees
% Utilisation du taux d'entropie

Ke=zeros(2);
Ke(1)=expand(u1e'*A*u1e);
Ke(4)=expand(u2e'*A*u2e);
Ke(2)=expand(u1e'*A*u2e);
Ke(3)=Ke(2);
Ke=Ke/V;

Kp=zeros(2);
Kp(1)=expand(u1p'*A*u1p);
Kp(4)=expand(u2p'*A*u2p);
Kp(2)=expand(u1p'*A*u2p);
Kp(3)=Kp(2);
Kp=Kp/V;

Kn=zeros(2);
Kn(1)=expand(u1n'*A*u1n);
Kn(4)=expand(u2n'*A*u2n);
Kn(2)=expand(u1n'*A*u2n);
Kn(3)=Kn(2);
Kn=Kn/V;
Kn=inv(Kn);
