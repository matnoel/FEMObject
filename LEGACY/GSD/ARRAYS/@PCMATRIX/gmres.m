function [u,varargout] = gmres(A,b,restart,tol,maxiter,varargin)
% function [u,varargout] = gmres(A,b,restart,tol,maxiter)
% resolution de Ax=b par gradient conjugue preconditionne
% b PCMATRIX ou PCRADIALMATRIX
% A PCMATRIX
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
if nargin<3
    restart=[];
end
if nargin<4
    tol=[];
end
if nargin<5
    maxiter = size(A,1);
end


if ~israndom(A)
    error('A deterministe : utiliser une autre methode de resolution')
end

[A,b,PC] = pcsystemupdate(A,b,varargin{:});

n=size(A,1);

b=double(PCMATRIX(b));


Amean = mean(A);
masse=getmasse(A);

precond = getcharin('precond',varargin,1);
switch precond
    case 1
        [ML,MU] = lu(Amean);
    case 2
        ML=cell(1,length(PC));
        MU=cell(1,length(PC));
        for p=1:length(ML)
            ML{p}=sparse(n,n);
            for k=1:length(masse)
                ML{p}= ML{p} + A.MULTIMATRIX{k}*masse{k}(p,p);
            end
            [ML{p},MU{p}]=lu(ML{p});
        end
end

u=zeros(size(b));
varargout=cell(1,nargout-1);
[u(:),varargout{:}]=gmres(@(X) Afun(X,A,n,length(PC)),b(:),restart,tol,maxiter,@(X) Mfun(X,ML,MU,n,length(PC),precond));

if nargout>1 && varargout{1}~=0
    warning('gmres not converged')
end
u=PCMATRIX(u,[n,1],PC);

function Y = Afun(X,A,n,P)

X=reshape(X,n,P);
Y=zeros(n,P);
masse=getmasse(A);

Y = multisum(A.MULTIMATRIX*(sparse(X)*masse));

Y=Y(:);

return

function Y = Mfun(X,ML,MU,n,P,precond)
X = reshape(X,n,P);
Y = zeros(n,P);
switch precond
    case 1
        Y=MU\(ML\X);
    case 2
        for k=1:length(ML)
            Y(:,k) = MU{k}\(ML{k}\X(:,k));
        end
end
Y = Y(:);
return

