function [u,result]=dsolve_radial_arnoldi(A,B,b,t,varargin)

fprintf('RESOLUTION PAR DECOMPOSITION SPECTRALE : ARNOLDI ... \n')
uref = getcharin('reference',varargin);
nbfoncmax = getcharin('nbfoncmax',varargin,10);
nbfoncmaxsimul = getcharin('nbfoncmaxsimul',varargin,10);
tol = getcharin('tol',varargin,1e-4);
orthocrit = getcharin('orthocrit',varargin,1e-8);
restart = getcharin('restart',varargin,0);

tolpcg = min(tol/100,1e-10);

infoiter = 'n';

if isa(A,'PCRADIAL') & isempty(A.Lmasse{1})
A=calc_ximasse(A);
end
if isa(B,'PCRADIAL') & isempty(B.Lmasse{1})
B=calc_ximasse(B);
end

PC = getPC(A);
P = getP(PC);
n=size(A,1);
nt = length(t)-1;
dt=t(2)-t(1);
global Mt
global Dt
Mt = speye(nt)*dt;
Dt = spdiags(ones(nt,1),0,nt,nt) + spdiags(-ones(nt-1,1),-1,nt,nt);
Dt = Mt\Dt;

result.error=zeros(1,nbfoncmax);
result.nfonctions = 0 ;

global Am
Am = mean(A);
Am = 1/2*(Am+Am');
Am = speye(size(Am,1));
u = PCRADIAL(PC);

bu = b;
bu = funV(bu,@(x)reshape(x,n,nt));

Vu = zeros(n,nt,0);

l0 = PCARRAY(zeros(1,P+1),PC);
l = PCARRAY(zeros(1,P+1),PC);
U = zeros(n,nt);
DU = zeros(n,nt);
U0 = zeros(n,nt);




for k=0:restart
if k>=1 
    fprintf('Actualisation #%d \n',k)
end
V = zeros(n,nt,0);
DV = zeros(n,nt,0);
l0(:,:)=1;l0(1,1)=1;
alpha0=norm(l0);l0=l0/alpha0;

for i=1:nbfoncmaxsimul

f0=expect(bu,l0);
A0 = expect(A,l0,l0);
B0 = expect(B,l0,l0);
[RL,RU] = lu(A0/dt+B0);
for j=1:nt
U(:,j)=RU\(RL\(f0(:,j)+(j>1)*A0/dt*U(:,j-(j>1))));
end


H0=sqrt(myprod(U,U));

for j=1:size(V,3)
U=U-V(:,:,j)*myprod(V(:,:,j),U);
end
H=sqrt(myprod(U,U));

fprintf('residu = %.3e',full(H/H0))

if H/H0<orthocrit
    disp(' -> break')
break
end

U = U / H ;
DU = U*Dt';
V(:,:,end+1)=U;
DV(:,:,end+1)=DU;

fprintf(' -> %d fonctions\n',size(V,3))


fU = funV(bu,@(x) myprod(U,x));
aU = funV(A*(DU),@(x) myprod(U,x));
bU = funV(B*U,@(x) myprod(U,x));
aU = aU+bU;

[l,flag] = gmres(aU,fU,[],tolpcg);

l0 = l ; U0 = U ;

end


fU = funV(bu,@(x) myprod(V,x));
ADV = funV(A,@(x) prod3(x,DV));
BV = funV(B,@(x) prod3(x,V));
aU = funV(ADV,@(x) myprod(V,x));
bU = funV(BV,@(x) myprod(V,x));
aU=aU+bU;

[l,flag] = gmres(aU,fU,[],tolpcg);

for j=1:size(V,3)
u = u + reshape(V(:,:,j),n*nt,1).*l(j,:);
bu = bu - (A*DV(:,:,j)+B*V(:,:,j))*l(j,:);
end


Vu(:,:,end+1:end+size(V,3)) = V; 
result.nfonctions = size(Vu,3);
if ischarin('reference',varargin)
result.error(i)=norm(u-expand(uref))/norm(expand(uref));
fprintf('%d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(i))
end


end

function v=myprod(U,V)

global Mt
for i=1:size(U,3)
for j=1:size(V,3)    
v(i,j)=sum(sum(U(:,:,i).*(V(:,:,j))));
end
end
return

function AV=prod3(A,V)

AV = zeros([size(A,1),size(V,2),size(V,3)]);
for j=1:size(V,3)    
AV(:,:,j)=A*V(:,:,j);
end

return

