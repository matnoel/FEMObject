%% Homogeneisation meca 2D separee en XY
% Definition des modeles
dim=2;
size_img=[1 1];
nelem=[50 60];

Dx=DOMAIN(1,0,size_img(1));
Dy=DOMAIN(1,0,size_img(2));

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
u=unfreevector(u1sep,1,Sxe);
u=unfreevector(u,2,Sye);
u=ubcx+u;

%% Export separe du deplacement

x=getcoord(getnode(Sxe));
y=getcoord(getnode(Sye));
nx=numel(x);
ny=numel(y);
zx=zeros(nx,1);
zy=zeros(ny,1);

nodes={x;y};

cx=(1:nx-1)';
cy=(1:ny-1)';
cells = {[cx cx+1];[cy cy+1]};

names={{'x'},{'y'}};

gu=gathervectors(u);
vector=gu.F{3};

ux=removedim(u,3);
ux.alpha=ux.alpha.*vector(1,:);
vx=normalizefuns(ux);
vx=normalizefactors(vx);
vx=gathervectors(vx)';

uy=removedim(u,3);
uy.alpha=uy.alpha.*vector(2,:);
vy=normalizefuns(uy);
vy=normalizefactors(vy);
vy=gathervectors(vy)';

nodes_fields = cell(2,2);
nodes_fields{1,1} = vx.F{1};
nodes_fields{2,1} = vx.F{2};

nodes_fields{1,2} = vy.F{1};
nodes_fields{2,2} = vy.F{2};

nodes_fields_names= cell(2,1);
nodes_fields_names{1}= 'ux';
nodes_fields_names{2}= 'uy';

cell_fields = cell(2,1);
cell_fields_name=cell(1,1);

filename= 'export_separe_dispXY.pxdmf';
output_mesh_PXDMF(filename, nodes, cells, names, nodes_fields, cell_fields,...
    nodes_fields_names, cell_fields_name, 'from1');

