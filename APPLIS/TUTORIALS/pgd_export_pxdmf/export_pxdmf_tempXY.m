%% Laplacien 2D separee en XY
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
[X,Y]=meshgrid( linspace(0,size_img(1),nelem(1)+1),...
    linspace(0,size_img(2),nelem(2)+1));
cx=0.5;
cy=0.5;
rfibre=0.3338;
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
[usep,resultref] = solve(Af,b,PGD);
ue=unfreevector(usep,1,Sxe);
ue=unfreevector(ue,2,Sye);
ue=ubcx+ue;

%% Export separe de la temperature

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

ve=normalizefuns(ue);
ve=normalizefactors(ve);
ve=gathervectors(ve)';

nodes_fields = cell(2,1);
nodes_fields{1,1} = ve.F{1};
nodes_fields{2,1} = ve.F{2};

nodes_fields_names= cell(1);
nodes_fields_names{1}= 'T';

cell_fields = cell(2,1);
cell_fields_name=cell(1,1);

filename= 'export_separe_tempXY.pxdmf';
output_mesh_PXDMF(filename, nodes, cells, names, nodes_fields, cell_fields,...
    nodes_fields_names, cell_fields_name, 'from1');
