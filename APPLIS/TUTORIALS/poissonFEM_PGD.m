% Some small examples to illustrate the use of FEMOBJECT
% 2D Laplacian problems in both regular and separated form
% 26/10/2010 G.L.



%% Small example, poisson, laplace equations
clear();

%Domain (geometry)
D=DOMAIN(2,[0,0],[1,1]);

%Mesh
nelem=30;
S=mesh(D,nelem,nelem);

%Create dofs
S=createddlnode(S,DDL('u'));



%Define Dirichlet BC
% S=addcl(S,getedge(D,1),'u',0);% u=0 on edge 1
% S=addcl(S,getedge(D,3),'u',1);% u=1 on edge 3
S=addcl(S,[],'u',0);% All the Bnd to Zero

% S=addcl(S,getedge(D,3),'u',@(x) x(:,1));% Non uniform u prescribed on edge 3

%Lin form
s=0;
x = getcoord(getnode(S));
s = 1+1.8*cos(pi*x(:,1)).*cos(2*pi*x(:,2));


l=LINFORM(0,s,0); %(Non)Uniform heat source inside the domain
% W : last zero is mandatory for the non-homogeneous case




%Bilin form
a=BILINFORM(1,1,1);%Bilinear form

%Build operators
A=a{S}(:,:);
L=l{S}(:);

%Solve system
sol=solve(A,L);

%Plot solution
plot(sol,S)
title('Solution, classical F.E.')

%% The same example, but solved using a XY PGD

clear();

% 1D Domain (geometry)
D1D=DOMAIN(1,0,1);

%Mesh
nelem=30;
S1D=mesh(D1D,nelem);

%Create dofs
S1D=createddlnode(S1D,DDL('u'));

%Define Dirichlet BC
%Case where u=0 on the whole bnd
S1D=addcl(S1D,POINT([0;1]),'u',0);% All the Bnd to Zero


%LinForm
s=1.;
l1D=LINFORM(0,s); %Uniform part of the heat source inside the domain

x = getcoord(getnode(S1D));
sx = 1.8*cos(pi*x(:,1));
sy = cos(2*pi*x(:,1));
l1Dx=LINFORM(0,sx,0); %NON Uniform part of the heat source inside the domain
l1Dy=LINFORM(0,sy,0); %NON Uniform part of the heat source inside the domain


%Bilinear Forms
k1D = BILINFORM(1,1);
m1D = BILINFORM(0,0);

%Assembly
K1D = k1D{S1D}(:,:);
M1D = m1D{S1D}(:,:);
L1D = l1D{S1D}{:};
L1Dx = l1Dx{S1D}{:};
L1Dy = l1Dy{S1D}{:};

Asep = SEPMATRIX({K1D,M1D;M1D,K1D});
bsep = SEPMATRIX({L1D,L1D}) + SEPMATRIX({L1Dx,L1Dy}) ;

%Solve
% Alternate solver, no update, metric given by the operator (Asep)
solver = SEPSOLVER(getdim(Asep),'residual',0,'maxorder',40,'tol',1e-5,'update',0,'metric',Asep);
[usep,result] = solve(Asep,bsep,solver);

figure(5)
semilogy(result.error,'g')
title('Convergence of the error w/r to the given metric')
xlabel('Nb modes')
ylabel('Error')
hold on

figure(6)
clf
temp = expand(usep);
[S,X,Y ] = mesh(DOMAIN(2),nelem,nelem);
S = createddlnode(S,DDL('u'));
S = addcl(S,[],'u');
plot(temp(:),S)%,'surface')
title('Solution, separated form')

%% Poisson with a hole...
clear();

%Domain (geometry)
D=DOMAIN(2,[0,0],[1,1]);

%Mesh
nelem=100;
S=mesh(D,nelem,nelem);

%Create Levelset
ls1 = LSCIRCLE(0.5,0.5,0.2);ls2 = LSCIRCLE(0.7,0.7,0.1);
ls= lseval(union(ls1,ls2),S);

figure(1)
clf
plot(ls,S,'surface')
contourplot(ls,S)
title('LevelSet + iso-zero')

%Create the ind. function
lsv = getvalue(ls); %get the vector of values of the lset
I = max(0,double(lsv>=0)); %Indic. func

figure(2)
clf
plot(I(:),S,'surface')
title('Indic function')

%Create dofs
S=createddlnode(S,DDL('u'));



%Define Dirichlet BC
% S=addcl(S,getedge(D,1),'u',0);% u=0 on edge 1
% S=addcl(S,getedge(D,3),'u',1);% u=1 on edge 3
S=addcl(S,[],'u',0);% All the Bnd to Zero

% S=addcl(S,getedge(D,3),'u',@(x) x(:,1));% Non uniform u prescribed on edge 3

%Lin form
s=0;
x = getcoord(getnode(S));
s = 1+1.8*cos(pi*x(:,1)).*cos(2*pi*x(:,2));


l=LINFORM(0,s.*I,0); %(Non)Uniform heat source inside the domain
% W : last zero is mandatory for the non-homogeneous case
% W : NOTE THE USE OF I !





%Bilin form
a=BILINFORM(1,1,I,0);%Bilinear form

%Build operators
A=a{S}(:,:);
L=l{S}(:);

%Solve system
sol=solve(A,L);

%Plot solution
figure(3)
clf
plot(sol,S)
title('Solution, classical F.E. with ind. function')

%% Poisson with a hole, separated form...

clear();

% 1D&2D Domains (geometry)
D1D=DOMAIN(1,0,1);
D2D=DOMAIN(2,[0,0],[1,1]); %Only needed to plot the LevelSet and the final solution

%Mesh
nelem=30;
S1D=mesh(D1D,nelem);
S2D=mesh(D2D,nelem,nelem);


%Create Levelset
ls1 = LSCIRCLE(0.5,0.5,0.2);ls2 = LSCIRCLE(0.7,0.7,0.1);
ls= lseval(union(ls1,ls2),S2D);
figure(1)
clf
plot(ls,S2D,'surface')
contourplot(ls,S2D)
title('LevelSet + iso-zero')

%Create ind. function
lsv = getvalue(ls); %get the vector of values of the lset
I = max(0,double(lsv>=0)); %Indic. func It should be smoothed for performances purpose.
I = reshape(I,nelem+1,nelem+1); %Mandatory to store the values in a nnodes x nnodes 2D xy array

% figure(2)
% clf
% plot(I(:),S2D,'surface')

%Separate I
[U,Sig,V]=svdtruncate(I,1e-4);
Im = U*Sig*V';%Approximated function with size(Sig) modes
figure(2)
clf
plot(Im(:),S2D,'surface')
plot(S2D)
title('Indic function')

%Create dofs (1D model)
S1D=createddlnode(S1D,DDL('u'));

%Define Dirichlet BC
%Case where u=0 on the whole bnd
S1D=addcl(S1D,POINT([0;1]),'u',0);% All the Bnd to Zero

%Because I has been replaced by a finite sum of x and  y modes :
%First, define the shape of Asep and bsep:
bsep = SEPMATRIX(2);
Asep = SEPMATRIX(2);

%Analyticaly separated source term
x = getcoord(getnode(S1D));
sx = 1.8*cos(pi*x(:,1));
sy = cos(2*pi*x(:,1));

%Loop on the modes of Im :
for i=1:size(Sig,1)

%LinForm (Uniform heat source part)
    s=1.;
    l1D_x=LINFORM(0,s*U(:,i),0); %To take into account x-mode N°i
    l1D_y=LINFORM(0,s*V(:,i),0); %To take into account y-mode N°i

    L1D_x = calc_vector(l1D_x,S1D); %Assemble the force vector
    L1D_y = calc_vector(l1D_y,S1D);

%LinForm (NON Uniform heat source part)
    l1Dx=LINFORM(0,sx.*U(:,i),0); %To take into account mode N°i
    l1Dy=LINFORM(0,sy.*V(:,i),0);

    L1Dx = calc_vector(l1Dx,S1D); %Assemble the force vector
    L1Dy = calc_vector(l1Dy,S1D);


%Bilinear Forms
    k1D = BILINFORM(1,1,U(:,i),0); %To take into account mode N°i
    m1D = BILINFORM(0,0,V(:,i),0);

    K1D = calc_matrix(k1D,S1D); %Assemble the stiffness matrix
    M1D = calc_matrix(m1D,S1D); %Assemble the mass matrix

%Accumulate the operators in Asep.
    Asep = Asep + SEPMATRIX({K1D,M1D;M1D,K1D},[Sig(i,i),Sig(i,i)]);% Scale it by singular values
    bsep = bsep + SEPMATRIX({L1D_x,L1D_y},Sig(i,i)) + SEPMATRIX({L1Dx,L1Dy},Sig(i,i)); % Scale it by singular values

end


%Solve
% Alternate solver, no update, metric given by the operator (Asep)
solver = SEPSOLVER(getdim(Asep),'residual',0,'maxorder',40,'tol',1e-5,'update',0,'metric',Asep);
[usep,result] = solve(Asep,bsep,solver);

%Convergence of the error
figure(5)
semilogy(result.error,'g')
title('Convergence of the error w/r to the given metric')
xlabel('Nb modes')
ylabel('Error')
hold on

figure(6)
clf

%Get the 2D solution :

% temp = expand(usep);

usep_unfree=usep;
usep_unfree=unfreevector(usep_unfree,1,S1D);
usep_unfree=unfreevector(usep_unfree,2,S1D);
temp = expand(usep_unfree);

S2D = createddlnode(S2D,DDL('u'));
S2D = addcl(S2D,[],'u');

plot(temp(:).*Im(:),S2D)
% plot(temp(:).*Im(:),S2D,'surface')
title('Solution, separated form + indic. function')

