%% Homogeneisation thermique 2D separee en XY, materiaux aleatoires
% Definition des modeles
dim=2;
size_img=[1 1];
r=30;
nelem=[r r+2];

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
indicf1=indicf1';

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
indicfs=permutedim(indicfs,[2;1]);

%% Introduction des variables aleatoires
km=1;
kf=10;

%rvm=RVNORMAL(km,0.1);
rvm=RVUNIFORM(km-0.1,km+0.3);
rvf=RVNORMAL(kf,2);
R=RANDVARS(rvm,rvf);

[X,PC]=PCTPMODEL(R,'order',5,'pcg');

oPC=one(PC);
oS=SEPMATRIX(oPC);
o=getphi(oPC);
o1=full(o{1}); % petits pbs avec le format sparse
o2=full(o{2});

B1 = calc_ximasse(X{1});
B1 = get_ximasse(B1);
B2 = calc_ximasse(X{2});
B2 = get_ximasse(B2);

%% Construction de l'operateur
kx=setfree(BILINFORM(1,1,1),0);
kx=kx{Sx}(:,:);
ky=setfree(BILINFORM(1,1,1),0);
ky=ky{Sy}(:,:);
mx=setfree(BILINFORM(0,0,1),0);
mx=mx{Sx}(:,:);
my=setfree(BILINFORM(0,0,1),0);
my=my{Sy}(:,:);
A1=SEPMATRIX({kx,my,B1{1},B1{2};...
    mx,ky,B1{1},B1{2}});

A2=SEPMATRIX(4);
A3=SEPMATRIX(4);
for i=1:indicfs.m
    kx=setfree(BILINFORM(1,1,indicfs.F{i,1},0),0);
    kx=kx{Sx}(:,:);
    ky=setfree(BILINFORM(1,1,indicfs.F{i,2},0),0);
    ky=ky{Sy}(:,:);
    mx=setfree(BILINFORM(0,0,indicfs.F{i,1},0),0);
    mx=mx{Sx}(:,:);
    my=setfree(BILINFORM(0,0,indicfs.F{i,2},0),0);
    my=my{Sy}(:,:);
    A2=A2+...
        SEPMATRIX({kx,my,B2{1},B2{2};...
        mx,ky,B2{1},B2{2}},...
        [indicfs.alpha(i),indicfs.alpha(i)]);
    A3=A3+...
        SEPMATRIX({kx,my,B1{1},B1{2};...
        mx,ky,B1{1},B1{2}},...
        [indicfs.alpha(i),indicfs.alpha(i)]);
end
A=A1+A2-A3;

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
ubcx=SEPMATRIX({x(:),ones(length(y),1),o1,o2}); 

b=-A*ubcx;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);

PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,...
    'update',1,'updatedim',3:4);
[u1sep,result1ref] = solve(Af,b,PGD);

u1e=unfreevector(u1sep,1,Sxe);
u1e=unfreevector(u1e,2,Sye);
u1e=ubcx+u1e;
u1e=PCTPMATRIX(u1e,PC,3:4);

% Mode 2
ubcy = SEPMATRIX({ones(length(x),1),y(:),o1,o2});
b=-A*ubcy;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,...
    'maxiter',10);
[u2sep,result2ref] = solve(Af,b,PGD);
u2e=unfreevector(u2sep,1,Sxe);
u2e=unfreevector(u2e,2,Sye);
u2e=ubcy+u2e;
u2e=PCTPMATRIX(u2e,PC,3:4);

%% Conditions aux limites periodiques
Sxp=Sx;
Syp=Sy;

Sxp=addclperiodic(Sxp,POINT(0),POINT(size_img(1)),'T');
Syp=addclperiodic(Syp,POINT(0),POINT(size_img(2)),'T');

lx = setfree( LINFORM(0,1),0);
Lx = calc_vector(lx,Sxp);
Ly = calc_vector(lx,Syp);
alpha=1;
Ap = A + SEPMATRIX({Lx*Lx',Ly*Ly',oS.F{2},oS.F{3}},alpha);

Af=freematrix(Ap,1,Sxp);
Af=freematrix(Af,2,Syp);

% Mode 1
b=-A*ubcx;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,'display',true,...
    'maxiter',10);
[u3sep,result3ref] = solve(Af,b,PGD);
u1p=unfreevector(u3sep,1,Sxp);
u1p=unfreevector(u1p,2,Syp);
u1p=ubcx+u1p;
u1p=PCTPMATRIX(u1p,PC,3:4);


% Mode 2
b=-A*ubcy;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,'display',true,...
    'maxiter',10);
[u4sep,result4ref] = solve(Af,b,PGD);
u2p=unfreevector(u4sep,1,Sxp);
u2p=unfreevector(u2p,2,Syp);
u2p=ubcy+u2p;
u2p=PCTPMATRIX(u2p,PC,3:4);

%% Conditions aux limites de Neumann
Sxn=Sx;
Syn=Sy;

l = setfree( LINFORM(0,1),0);
Lx = calc_vector(l,Sxn);
Ly = calc_vector(l,Syn);
alpha=1;
Ap = A + SEPMATRIX({Lx*Lx',Ly*Ly',oS.F{2},oS.F{3}},alpha);

% Mode 1
xm=zeros(nelem(1)+1,1);
xm(1)=-1;
xp=zeros(nelem(1)+1,1);
xp(nelem(1)+1)=1;
b = SEPMATRIX({xm,Ly,o1,o2;xp,Ly,o1,o2});
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,'update',1,...
    'updatedim',3:4,'maxiter',10);
[u5sep,result5ref] = solve(Ap,b,PGD);
u1n=u5sep;
u1n=PCTPMATRIX(u1n,PC,3:4);

% Mode 2
ym=zeros(nelem(2)+1,1);
ym(1)=-1;
yp=zeros(nelem(2)+1,1);
yp(nelem(2)+1)=1;
b = SEPMATRIX({Lx,ym,o1,o2;Lx,yp,o1,o2});
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-3,'display',true,...
    'maxiter',10);
[u6sep,result6ref] = solve(Ap,b,PGD);
u2n=u6sep;
u2n=PCTPMATRIX(u2n,PC,3:4);

%% Determination des conductivites homogeneisees
% Integration du gradient de temperature et du flux sur le domaine 2d

nbreal=50;
x=zeros(nbreal,2);
kmat=zeros(nbreal,1);
kfib=zeros(nbreal,1);

Ke=zeros(nbreal,4);
Kp=zeros(nbreal,4);
Kn=zeros(nbreal,4);

kk=zeros(2);

L=setfree(LINFORM(1,1),0);
L=L{S}(:);

for i=1:nbreal
    [kmat(i),x(i,:)]=random(X{1});
    kfib(i)=randomeval(X{2},x(i,:));

    kmap=kmat(i)+(kfib(i)-kmat(i))*indicf1';
    kmap=kmap(:);

    Lq=setfree(LINFORM(1,kmap,0),0);
    Lq=Lq{S}(:);

    u1=full(randomeval(u1e,x(i,:)))';
    u2=full(randomeval(u2e,x(i,:)))';    
    u1=u1(:);
    u2=u2(:);
    Ke(i,1:2)=Lq'*u1;
    Ke(i,3:4)=Lq'*u2;

    u1=randomeval(u1p,x(i,:))';
    u2=randomeval(u2p,x(i,:))';
    u1=u1(:);
    u2=u2(:);
    Kp(i,1:2)=Lq'*u1;
    Kp(i,3:4)=Lq'*u2;

    u1=randomeval(u1n,x(i,:))';
    u2=randomeval(u2n,x(i,:))';
    u1=u1(:);
    u2=u2(:);
    kk(1:2)=L'*u1;
    kk(3:4)=L'*u2;
    kk=inv(kk);
    Kn(i,1:2)=kk(1:2);
    Kn(i,3:4)=kk(3:4);
end

Ke=Ke/V;
Kp=Kp/V;
Kn=Kn/V;

K=zeros(nbreal,3);
K(:,1)=Ke(:,1);
K(:,2)=Kp(:,1);
K(:,3)=Kn(:,1);

disp([kmat kfib K])
