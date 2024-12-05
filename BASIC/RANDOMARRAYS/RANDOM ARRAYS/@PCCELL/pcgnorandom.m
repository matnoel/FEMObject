function [u,varargout] = pcgnorandom(A,b,tol,maxiter,varargin)
% function [u,varargout] = pcgnorandom(A,b,tol,maxiter)
% resolution de Ax=b par gradient conjugue preconditionne
% b PCCELL ou PCARRAY ou PCRADIAL
% tol : precision souhaitee (1e-13 par defaut ou si tol=[])
% maxiter : maximum iteration (100 par defaut ou si maxiter=[])
%
% function [u,varargout] = pcg(A,b,tol,maxiter,PC)
% PC : POLYCHAOS sur lequel est defini b et la solution x (A peut etre defini sur un autre)
%
% function [u,varargout] = pcg(A,b,tol,maxiter,'norandom',norandom)
% les composantes norandom de x sont deterministes
%
% function [u,varargout] = pcg(A,b,tol,maxiter,'precond',precond)
% precond=1 : preconditionneur moyenne (par defaut)
% precond=2 : preconditionneur bloc

n=size(A,1);
nosto = getcharin('norandom',varargin);
nns = length(nosto);
sto = setdiff(1:n,nosto);
ns = length(sto);


if nargin<3 | isempty(tol)
    tol=1e-13;
end
if nargin<4 | isempty(maxiter)
    maxiter = 100;
end

if nns==0
    varargout=cell(1,nargout-1);
    [u,varargout{:}] = pcg(A,b,tol,maxiter,varargin{:});
    return
end

PC = getclassin('POLYCHAOS',varargin,getPC(A));
PC = calc_masse(PC,getPC(A));
P = getP(PC);

if isa(b,'double') 
b=b.*one(PC);
end
b=double(expand(b));

b1 = b(sto,:);
b2 = b(nosto,:)*double(one(PC))';
b=[b1(:);b2(:)];


global masse
global masseun
masse=getmasse(PC);
for k=1:length(masse)
    masseun{k} = masse{k}*double(one(PC))';
end
global Amean;
Amean = mean(A);
MLnosto = chol(Amean(nosto,nosto));

precond = getcharin('precond',varargin,1);
switch precond
    case 1
MLsto = chol(Amean(sto,sto));
    case 2
MLsto=cell(1,length(PC));
for p=1:length(MLsto)
MLsto{p}=sparse(ns,ns);
for k=1:length(masse)
MLsto{p}= MLsto{p} + A.value{k}(sto,sto)*masse{k}(p,p);
end
MLsto{p}=chol(MLsto{p});
end
end

varargout=cell(1,nargout-1);

[v,varargout{:}]=pcg(@(X) Afun(X,A,length(PC),n,ns,nns,sto,nosto),...
      b(:),tol,maxiter,@(X) Mfun(X,MLsto,MLnosto,length(PC),n,ns,nns,sto,nosto,precond));
u=zeros(n,P+1);
v1 = v(1:ns*(P+1));
v2 = v(ns*(P+1)+1:end);

u(sto,:)=reshape(v1,ns,P+1);
u(nosto,:)=v2*double(one(PC));

u=PCARRAY(u,PC);

function Y = Afun(X,A,P,n,ns,nns,sto,nosto)

X1=reshape(X(1:ns*P),ns,P);
X2=reshape(X(ns*P+1:end),nns,1);
Y1=zeros(ns,P);
Y2=zeros(nns,1);
global masse
global masseun
global Amean
for k=1:length(masse)
Y1=Y1+A.value{k}(sto,sto)*X1*masse{k};
Y1=Y1+A.value{k}(sto,nosto)*X2*masseun{k}';
Y2 = Y2 + A.value{k}(nosto,sto)*X1*masseun{k};
end
Y2 = Y2 + Amean(nosto,nosto)*X2;
Y=[Y1(:);Y2(:)];
return

function Y = Mfun(X,MLsto,MLnosto,P,n,ns,nns,sto,nosto,precond)

X1=reshape(X(1:ns*P),ns,P);
X2=reshape(X(ns*P+1:end),nns,1);
Y1=zeros(ns,P);
Y2=zeros(nns,1);
switch precond
    case 1
        for k=1:P
Y1 = MLsto\(MLsto'\X1);
        end
    case 2
        for k=1:P
Y1(:,k) = MLsto{k}\(MLsto{k}'\X1(:,k));
        end        
end

Y2 = MLnosto\(MLnosto'\X2);

Y=[Y1(:);Y2(:)];
return

