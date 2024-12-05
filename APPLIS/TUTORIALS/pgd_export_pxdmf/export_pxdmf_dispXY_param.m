%% Mécanique, separation XY parametrique sur la force volumique
% Carre encastre sur ses bords
% Une force volumique est imposee parametre par un angle t

% Construction du modele

tmin=0;
tmax=2*pi;

dim=2;
size_img=[1 1];
nelem=[100 100];

Dx=DOMAIN(1,0,size_img(1));
Dy=DOMAIN(1,0,size_img(2));
Dt=DOMAIN(1,tmin,tmax);

Sx=mesh(Dx,nelem(1));
Sy=mesh(Dy,nelem(2));
St=mesh(Dt,180);

Sx=createddlnode(Sx,DDL('u'));
Sy=createddlnode(Sy,DDL('u'));
St=createddlnode(St,DDL('u'));

Sx=addcl(Sx,POINT([0;1]),'u');
Sy=addcl(Sy,POINT([0;1]),'u');

E=1;
nu=0.3;
C=elasticity_tensor(E,nu,2);

%% Construction de l'operateur

Ex=zeros(3,2);
Ex(1,1)=1;
Ex(3,2)=1;

Ey=zeros(3,2);
Ey(2,2)=1;
Ey(3,1)=1;

kx=BILINFORM(1,1,1);
kx=kx{Sx}(:,:);
mx=BILINFORM(0,0,1);
mx=mx{Sx}(:,:);
dxl=BILINFORM(1,0,1);
dxl=dxl{Sx}(:,:);
dxr=dxl';

ky=BILINFORM(1,1,1);
ky=ky{Sy}(:,:);
my=BILINFORM(0,0,1);
my=my{Sy}(:,:);
dyl=BILINFORM(1,0,1);
dyl=dyl{Sy}(:,:);
dyr=dyl';

mt=BILINFORM(0,0,1);
mt=mt{St}(:,:);

Gxx=Ex'*C*Ex;
Gyy=Ey'*C*Ey;
Gxy=Ex'*C*Ey;
Gyx=Gxy';

A=SEPMATRIX({kx,my,Gxx,mt;mx,ky,Gyy,mt;dxl,dyr,Gxy,mt;dxr,dyl,Gyx,mt});

%% Construction du 2nd membre

x=getcoord(getnode(Sx));
y=getcoord(getnode(Sy));
t=getcoord(getnode(St));

wx=exp(-50*(x-0.5).^2);
wy=exp(-50*(y-0.5).^2);
g1=[1;0];
g2=[0;1];
wt1=cos(t);
wt2=sin(t);

lx=LINFORM(0,wx,0);
lx=lx{Sx}(:);
ly=LINFORM(0,wy,0);
ly=ly{Sy}(:);
lt1=LINFORM(0,wt1,0);
lt1=lt1{St}(:);
lt2=LINFORM(0,wt2,0);
lt2=lt2{St}(:);

b1=SEPMATRIX({lx,ly,g1,lt1});
b2=SEPMATRIX({lx,ly,g2,lt2});

b=b1+b2;

%% Resolution

PGD = SEPSOLVER(getdim(A),'maxorder',100,'tol',1e-5);
[u,result] = solve(A,b,PGD);

u=full(u);

u=unfreevector(u,1,Sx);
u=unfreevector(u,2,Sy);
u=unfreevector(u,4,St);

%% Export separe parametrique

nx=numel(x);
ny=numel(y);
nt=numel(t);

nodes={x;y;t};

cx=(1:nx-1)';
cy=(1:ny-1)';
ct=(1:nt-1)';
cells = {[cx cx+1];[cy cy+1];[ct ct+1]};

names={{'x'},{'y'},{'t'}};

un=normalizefuns(u);
un=normalizefactors(un);

gu=gathervectors(un);
vector=gu.F{3};

vx=removedim(un,3);
vx=removedim(vx,3);
vx.alpha=vx.alpha.*vector(1,:);
vx=normalizefactors(vx);
vx=gathervectors(vx)';

vy=removedim(un,3);
vy=removedim(vy,3);
vy.alpha=vy.alpha.*vector(2,:);
vy=normalizefactors(vy);
vy=gathervectors(vy)';

vt=gu.F{4}';

nodes_fields = cell(3,2);
nodes_fields{1,1} = vx.F{1};
nodes_fields{2,1} = vx.F{2};

nodes_fields{1,2} = vy.F{1};
nodes_fields{2,2} = vy.F{2};

nodes_fields{3,1} = vt;
nodes_fields{3,2} = vt;

nodes_fields_names= cell(3,1);
nodes_fields_names{1}= 'ux';
nodes_fields_names{2}= 'uy';
nodes_fields_names{3}= 't';

cell_fields = cell(3,1);
cell_fields_name=cell(1,1);

filename= 'export_separe_dispXY_param.pxdmf';
output_mesh_PXDMF(filename, nodes, cells, names, nodes_fields, cell_fields,...
    nodes_fields_names, cell_fields_name, 'from1');




