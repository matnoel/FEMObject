
rx=40;
rt=50;
dims = 2;
calcref = 1;

T = TIMEMODEL(0,.03,rt);
N = DGTIMESOLVER(T,0);
Mt = getMmatrix(N);
Dt = getDmatrix(N,'basic');
nt = length(N);
ft = getMmatrix(N)*full(double(one(N))');

Sx = mesh(DOMAIN(1,0,0.7),ceil(rx*.7));
Sx = union(Sx,mesh(DOMAIN(1,0.7,0.8),ceil(rx*.1)));
Sx = union(Sx,mesh(DOMAIN(1,0.8,1),ceil(rx*.2)));
Sx = createddlnode(Sx,DDL('u'),DDL('g'));
Sx = addcl(Sx,POINT([0;1]),'u',0);
nx = getnbddlfree(Sx);
xnode = getcoord(getnode(Sx)); 
ynode = getcoord(getnode(Sx)); 
vxnodex = 1;
vxnodey = -(ynode-1/2);
vynodex = -(1/2-xnode);
vynodey = 1;
ax = BILINFORM(1,1);
cxx = BILINFORM(0,1,vxnodex,0);
cxy = BILINFORM(0,0,vxnodey,0);
cyx = BILINFORM(0,0,vynodex,0);
cyy = BILINFORM(0,1,vynodey,0);
mx = BILINFORM(0,0);
Ax = ax{Sx}(:,:);
Ay = Ax;
Mx = mx{Sx}(:,:);
My = Mx;
Cxx = cxx{Sx}(:,:);
Cxy = cxy{Sx}(:,:);
Cyx = cyx{Sx}(:,:);
Cyy = cyy{Sx}(:,:);

lx = LINFORM(0,1);
fx = lx{keepgroupelem(Sx,2)}(:);
fy = fx;
Sy = Sx;
nx = getnbddlfree(Sx);
ny = getnbddlfree(Sy);


ka=1;
ca=250;
ra=10;


A = SEPMATRIX({Mx,My,Dt;ka*Ax,My,Mt;Mx,ka*Ay,Mt});
A = A + SEPMATRIX({ca*Cxx,Cxy,Mt;Cyx,ca*Cyy,Mt}) ...
    + SEPMATRIX({ra*Mx,My,Mt});
b = SEPMATRIX({fx,fy,ft});


%%

[S,X,Y] = mesh(DOMAIN(2),xnode,ynode);
S =  createddlnode(S,DDL('u'),DDL('g'));
S = addcl(S,[],'u',0);

vx = @(x) (x(:,2)-1/2);
vy = @(x) (1/2-x(:,1));
nodeco = getcoord(getnode(S)); 
VX = vx(nodeco);
VY = vy(nodeco);

a = BILINFORM(1,1);
c = BILINFORM(0,1,[VX,VY],0);
m = BILINFORM(0,0);
AS = ka*a{S}(:,:)+ca*c{S}(:,:) + ra*m{S}(:,:);
MS = m{S}(:,:);
l = LINFORM(0,1);
zz=intersect(S,DOMAIN(2,[0.7,0.7],[0.8,0.8]));
fS = l{zz}(:);


if calcref
    N = setparam(N,'display',true)
    uref = dsolve(N,fS*one(N),MS,AS);
    uref = unfreevector(S,uref);

    s=[size(X,1),size(X,2),length(N)];    
    solref = reshape(double(uref),s);
    solref = freevector(Sx,solref,1);
    solref = freevector(Sy,solref,2);
end


% affichage reference 2D
gg=reshape(solref,size(solref,1)*size(solref,2),size(solref,3));
figure(45)
clf
evol(TIMEMATRIX(gg,N),S)
%% SEPARE ESPACE-ESPACE-TEMPS
if calcref==0
    solref=[];
end
minres=0;
maxorder = 40;
maxiter=10;
update = 0;
alphaupdate = 1;
updatedim = 1:3;
adjoint = 0;
tol=1e-10
updatestep=1;
itercritupdate=1e-1;
updateadjoint=false;
adjointtype=0;
if adjoint && minres==0
    col = 'r';
elseif adjoint && minres
    col = 'y';    
elseif minres
    col = 'g';
else
    col = 'b';
end
if update
    col=[col,'-^'];
else
    col=[col,'-*'];
end

errorindicator = 'residual';



solver = SEPSOLVER(getdim(A),'tol',tol,'alphaupdate',alphaupdate,...
    'update',update,'updatedim',updatedim,...
    'maxorder',maxorder,'maxiter',maxiter,'reference',solref,...
    'adjoint',adjoint,'updatestep',updatestep,...
    'itercritupdate',itercritupdate,...
    'errorindicator',errorindicator,'updateadjoint',updateadjoint,...
    'adjointtype',adjointtype);
if minres==1
    [u,result] = solve(A'*A,A'*b,solver);
else
    [u,result] = solve(A,b,solver);    
end
u = unfreevector(unfreevector(u,2,Sy),1,Sx);
alpha = getalpha(u);
Vx = extractvectors(u,1);
Vy = extractvectors(u,2);
Vt = extractvectors(u,3);


figure(17);semilogy(result.error,col);
hold on

%%

%% AFFICHE 2D
Vxy = zeros(numel(X),length(alpha));

for i=1:length(alpha)
    vtemp = Vx(:,i)*Vy(:,i)';    
    Vxy(:,i)=vtemp(:);
end

ugsd = Vxy*diag(alpha)*Vt';
ugsd = TIMEMATRIX(ugsd,N,[size(ugsd,1),1]);


figure(5)
evol(ugsd,S,'rescale','z')

%% SEPARE ESPACE-TEMPS
Abis = SEPMATRIX({MS,Dt;AS,Mt});
bbis = SEPMATRIX({fS,ft});

if calcref==0
    solrefbis=[];
else
    solrefbis = freevector(S,uref);
end
maxorder = 20;
maxiter=10;
update = 0;
alphaupdate = 0;
updatedim = 2;
adjoint =0;
tol=1e-10
updatestep=1;
itercritupdate=1e-1;
if adjoint
    col = 'r';
else
    col = 'b';    
end
if update
    col=[col,'-^'];
else
    col=[col,'-*'];
end

errorindicator = 'residual';


solver = SEPSOLVER(getdim(Abis),'tol',tol,'alphaupdate',alphaupdate,...
    'update',update,'updatedim',updatedim,...
    'maxorder',maxorder,'maxiter',maxiter,'reference',solrefbis,...
    'adjoint',adjoint,'updatestep',updatestep,...
    'itercritupdate',itercritupdate,...
    'errorindicator',errorindicator,'updateadjoint',false);
[u,result] = solve(Abis,bbis,solver);

u = unfreevector(u,1,S);
alpha = getalpha(u);
Vx = extractvectors(u,1);
Vt = extractvectors(u,2);
ut = TIMEMATRIX(Vx*diag(alpha)*Vt',N);   


figure(17);semilogy(result.error,col);
hold on

%%
figure(65)
clf
evol(ut,S)
