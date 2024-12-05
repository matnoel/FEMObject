function [u,varargout] = pcg(A,b,tol, maxiter,varargin)

Amean = mean(A);
ML = chol(Amean);

if nargin<3 | isempty(tol)
    tol=1e-13;
end
if nargin<4 | isempty(maxiter)
    maxiter = 100;
end

A=calc_ximasse(A);

if isa(b,'double') 
b=b.*one(A);
end
b=double(expand(b));

n=size(b,1);
P = getP(A);

u=zeros(size(b));
varargout=cell(1,nargout-1);

[u(:),varargout{:}]=pcg(@(X) Afun(X,A,n,P),b(:),tol,maxiter,@(X) Mfun(X,ML,n,P));

u=PCARRAY(u,getPC(A));

function Y = Afun(X,A,n,P)

X=reshape(X,n,P+1);
Y=zeros(n,P+1);

for k=1:length(A)
Y=Y+A.V{k}*X*A.Lmasse{k};
end
Y=Y(:);

return

function Y = Mfun(X,ML,n,P)

Y=ML\(ML'\reshape(X,n,P+1));
Y = Y(:);

return

