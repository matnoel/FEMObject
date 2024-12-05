function [u,varargout] = gmres(A,b,restart,tol, maxiter,varargin)
if nargin<3
    restart=[];
end
if nargin<4
    tol=[];
end
if nargin<5
    maxiter = [];
end

Amean = mean(A);
[ML,MU] = lu(Amean);

if isempty(A.Lmasse{1}) | length(A.Lmasse)~=A.m
A=calc_ximasse(A);
end



if isa(b,'double') 
b=b.*one(A);
end
b=double(expand(b));


n=size(b,1);
P = getP(A);

u=zeros(size(b));
varargout=cell(1,nargout-1);

[u(:),varargout{:}]=gmres(@(X) Afun(X,A,n,P),b(:),restart,max([tol,1e-12]),max([maxiter,100]),@(X) Mfun(X,ML,MU,n,P));

u=PCARRAY(u,getPC(A));

function Y = Afun(X,A,n,P)

X=reshape(X,n,P+1);
Y=zeros(n,P+1);

for k=1:length(A)
Y=Y+A.V{k}*X*A.Lmasse{k};
end
Y=Y(:);

return

function Y = Mfun(X,ML,MU,n,P)

Y=MU\(ML\reshape(X,n,P+1));
Y = Y(:);   
return

