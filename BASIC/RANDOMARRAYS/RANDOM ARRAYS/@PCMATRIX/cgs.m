function [u,result] = cgs(A,b,tol,maxiter,varargin)
% function [u,varargout] = cgs(A,b,tol,maxiter)
% resolution de Ax=b par gradient conjugue preconditionne
% b PCMATRIX ou PCRADIALMATRIX
% A PCMATRIX
% tol : precision souhaitee (1e-13 par defaut ou si tol=[])
% maxiter : maximum iteration (100 par defaut ou si maxiter=[])
%
% function [u,varargout] = cgs(A,b,tol,maxiter,PC)
% PC : POLYCHAOS sur lequel est defini b et la solution x (A peut etre
% defini sur un autre)
%
% function [u,varargout] = cgs(A,b,tol,maxiter,'precond',precond)
% precond=1 : preconditionneur moyenne (par defaut)
% precond=2 : preconditionneur bloc


if nargin<3
    tol=[];
end
if nargin<4
    maxiter = 1000;
end


if ~israndom(A)
    error('A deterministe : utiliser une autre methode de resolution')
end

[A,b,PC] = pcsystemupdate(A,b,varargin{:});

clock0=clock;


n=size(A,1);

b=double(cell2mat(PCMATRIX(b)));


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

[u(:),flag,relres,iter,resvec]=cgs(@(X) Afun(X,A,n,length(PC)),b(:),tol,maxiter,@(X) Mfun(X,ML,MU,n,length(PC),precond));

u=PCMATRIX(u,[n,1],PC);

if flag~=0
    warning('cgs not converged')
end


if nargout>1
    result.totaltime = etime(clock,clock0);
    result.flag = flag;
    result.relres = relres;
    result.iter = iter;
    result.resvec = resvec;
end


function Y = Afun(X,A,n,P)

X=reshape(X,n,P);
Y=zeros(n,P);
masse=getmasse(A);

temp = sparse(X)*masse;

if iscell(A.MULTIMATRIX)
    temp = mat2cell(temp);
end
Y = multisum(A.MULTIMATRIX*temp);

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

