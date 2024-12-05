%% Homogeneisation meca 2D separee en XY
% Definition des modeles
dim=2;
size_img=[1 1];
nelem=[50 60];

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

%% Separation de l'indicatrice et creation de la carte de materiaux
Ef=10;
nuf=0.3;
Cf=elasticity_tensor(Ef,nuf,2);
Em=1;
num=0.3;
Cm=elasticity_tensor(Em,num,2);
Cd=Cf-Cm;

C1=SEPMATRIX(3);
C1.alpha=1;
C1.F{1,1}=ones(nelem(1)+1,1);
C1.F{1,2}=ones(nelem(2)+1,1);
C1.F{1,3}=Cm;

indicfs=multisvd(indicf,'maxorder',50,'tol',1e-3);
indicfs=permutedim(indicfs,[2;1]); % X et Y sont remis dans l'ordre
C2=SEPMATRIX(3);
C2.alpha=indicfs.alpha;
C2.F=indicfs.F;
for i=1:indicfs.m
    C2.F{i,3}=Cd;
end

Cmap=C1+C2;

%% Construction de l'operateur
Ex=zeros(3,2);
Ex(1,1)=1;
Ex(3,2)=1;

Ey=zeros(3,2);
Ey(2,2)=1;
Ey(3,1)=1;

A=SEPMATRIX(3);
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

    A=A+SEPMATRIX({kx,my,Gxx;mx,ky,Gyy;dxl,dyr,Gxy;dxr,dyl,Gyx},...
        [Cmap.alpha(i),Cmap.alpha(i),Cmap.alpha(i),Cmap.alpha(i)]);
end

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
ubcx = SEPMATRIX({x(:),ones(length(y),1),[1;0]});
b=-A*ubcx;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u1sep,result1ref] = solve(Af,b,PGD);
u1e=unfreevector(u1sep,1,Sxe);
u1e=unfreevector(u1e,2,Sye);
u1e=ubcx+u1e;

% Mode 2
ubcy = SEPMATRIX({ones(length(x),1),y(:),[0;1]});
b=-A*ubcy;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u2sep,result2ref] = solve(Af,b,PGD);
u2e=unfreevector(u2sep,1,Sxe);
u2e=unfreevector(u2e,2,Sye);
u2e=ubcy+u2e;

% Mode 3
ubcm = SEPMATRIX({x(:),ones(length(y),1),[0;1]});
b=-A*ubcm;
b=freevector(b,1,Sxe);
b=freevector(b,2,Sye);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u3sep,result3ref] = solve(Af,b,PGD);
u3e=unfreevector(u3sep,1,Sxe);
u3e=unfreevector(u3e,2,Sye);
u3e=ubcm+u3e;

%% Conditions aux limites periodiques
Sxp=Sx;
Syp=Sy;

Sxp=addclperiodic(Sxp,POINT(0),POINT(size_img(1)),'u');
Syp=addclperiodic(Syp,POINT(0),POINT(size_img(2)),'u');

lx = setfree( LINFORM(0,1),0);
Lx = calc_vector(lx,Sxp);
Ly = calc_vector(lx,Syp);
alpha=1;
Ap = A + SEPMATRIX({Lx*Lx',Ly*Ly',eye(2)},alpha);

Af=freematrix(Ap,1,Sxp);
Af=freematrix(Af,2,Syp);

% Mode 1
b=-A*ubcx;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u4sep,result4ref] = solve(Af,b,PGD);
u1p=unfreevector(u4sep,1,Sxp);
u1p=unfreevector(u1p,2,Syp);
u1p=ubcx+u1p;

% Mode 2
b=-A*ubcy;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u5sep,result5ref] = solve(Af,b,PGD);
u2p=unfreevector(u5sep,1,Sxp);
u2p=unfreevector(u2p,2,Syp);
u2p=ubcy+u2p;

% Mode 3
b=-A*ubcm;
b=freevector(b,1,Sxp);
b=freevector(b,2,Syp);
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u6sep,result6ref] = solve(Af,b,PGD);
u3p=unfreevector(u6sep,1,Sxp);
u3p=unfreevector(u3p,2,Syp);
u3p=ubcm+u3p;

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
b = SEPMATRIX({xm,ly,[1;0];xp,ly,[1;0]});
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u7sep,result7ref] = solve(Ap,b,PGD);
u1n=u7sep;

% Mode 2
lx=calc_vector(l,Sxn);
ym=zeros(nelem(2)+1,1);
ym(1)=-1;
yp=zeros(nelem(2)+1,1);
yp(nelem(2)+1)=1;
b = SEPMATRIX({lx,ym,[0;1];lx,yp,[0;1]});
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u8sep,result8ref] = solve(Ap,b,PGD);
u2n=u8sep;

% Mode 3
b = SEPMATRIX({xm,ly,[0;1];xp,ly,[0;1];...
    lx,ym,[1;0];lx,yp,[1;0]});
PGD = SEPSOLVER(getdim(A),'maxorder',30,'tol',1e-5);
[u9sep,result9ref] = solve(Ap,b,PGD);
u3n=u9sep;

%% Determination du tenseur d'elasticite homogeneise
% Utilisation de l'energie de deformation

u={u1e,u2e,u3e};

Ce=zeros(3);
for i=1:3
    for j=1:3
        Ce(i,j)=expand(u{i}'*A*u{j});
    end
end
Ce=Ce/V;

u={u1p,u2p,u3p};

Cp=zeros(3);
for i=1:3
    for j=1:3
        Cp(i,j)=expand(u{i}'*A*u{j});
    end
end
Cp=Cp/V;

u={u1n,u2n,u3n};

Cn=zeros(3);
for i=1:3
    for j=1:3
        Cn(i,j)=expand(u{i}'*A*u{j});
    end
end
Cn=Cn/V;
Cn=inv(Cn);
