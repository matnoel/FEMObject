function [u,varargout] = gmres(A,b,restart,tol, maxiter,varargin)

if nargin<3
    restart=[];
end
if nargin<4 
    tol=[];
end
if nargin<5
    maxiter = size(A,1);
end

varargout=cell(1,nargout-1);

if isa(A,'PCMATRIX')
 [u,varargout{:}]=gmres(A,b,restart,tol,maxiter,varargin{:});  
 return
end

if ~israndom(A)
    error('A deterministe : utiliser une autre methode de resolution')
end

Amean = mean(A);
[ML,MU] = lu(Amean);

[A,b,PC]=pcsystemupdate(A,b,varargin{:});

n=size(A,1);

b=double(PCMATRIX(b));

u=zeros(size(b));
varargout=cell(1,nargout-1);

[u(:),varargout{:}]=gmres(@(X) Afun(X,A,n,length(PC)),b(:),restart,tol,maxiter,@(X) Mfun(X,ML,MU,n,length(PC)));


u=PCMATRIX(u,[n,1],PC);

function Y = Afun(X,A,n,P)

X=reshape(X,n,P);
Y=zeros(n,P);

Y = A.V*(X*A.DLmasse);
Y = multisum(Y);
Y=Y(:);

return

function Y = Mfun(X,ML,MU,n,P)

Y=MU\(ML\reshape(X,n,P));
Y = Y(:);   
return

