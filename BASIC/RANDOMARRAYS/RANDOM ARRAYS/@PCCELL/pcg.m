function [u,varargout] = pcg(A,b,tol,maxiter,varargin)
% function [u,varargout] = pcg(A,b,tol,maxiter)
% resolution de Ax=b par gradient conjugue preconditionne
% b PCCELL ou PCARRAY ou PCRADIAL
% A PCCELL
% tol : precision souhaitee (1e-13 par defaut ou si tol=[])
% maxiter : maximum iteration (100 par defaut ou si maxiter=[])
%
% function [u,varargout] = pcg(A,b,tol,maxiter,PC)
% PC : POLYCHAOS sur lequel est defini b et la solution x (A peut etre
% defini sur un autre)
%
% function [u,varargout] = pcg(A,b,tol,maxiter,'precond',precond)
% precond=1 : preconditionneur moyenne (par defaut)
% precond=2 : preconditionneur bloc


PC = getclassin('POLYCHAOS',varargin,getPC(A));
PC = calc_masse(PC,getPC(A));
P = getP(PC);


if nargin<3 | isempty(tol)
    tol=1e-13;
end
if nargin<4 | isempty(maxiter)
    maxiter = 100;
end

if isa(b,'double') 
b=b.*one(PC);
end
b=double(expand(b));
n=size(b,1);

global masse
masse=getmasse(PC);
Amean = mean(A);
ML = chol(Amean);
precond = getcharin('precond',varargin,1);
switch precond
    case 1
ML = chol(Amean);
    case 2
ML=cell(1,length(PC));
for p=1:length(ML)
ML{p}=sparse(n,n);
for k=1:length(masse)
ML{p}= ML{p} + A.value{k}*masse{k}(p,p);
end
ML{p}=chol(ML{p});
end
end


u=zeros(size(b));
varargout=cell(1,nargout-1);

[u(:),varargout{:}]=pcg(@(X) Afun(X,A,n,length(PC)),b(:),tol,maxiter,@(X) Mfun(X,ML,n,length(PC),precond));

u=PCARRAY(u,PC);

function Y = Afun(X,A,n,P)

X=reshape(X,n,P);
Y=zeros(n,P);
global masse

for k=1:length(masse)
Y=Y+A.value{k}*X*masse{k};
end
Y=Y(:);
return

function Y = Mfun(X,ML,n,P,precond)
X = reshape(X,n,P);
Y = zeros(n,P);
switch precond
    case 1
        Y=ML\(ML'\X);
    case 2
   for k=1:length(ML)
   Y(:,k) = ML{k}\(ML{k}'\X(:,k));
   end
end
Y = Y(:);
return

