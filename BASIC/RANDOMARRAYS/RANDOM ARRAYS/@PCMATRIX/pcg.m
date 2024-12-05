function [u,result] = pcg(A,b,tol,maxiter,varargin)
% function [u,varargout] = pcgnorandom(A,b,tol,maxiter)
% resolution de Ax=b par gradient conjugue preconditionne
% b PCMATRIX ou PCRADIALMATRIX
% A PCMATRIX
% tol : precision souhaitee (1e-6 par defaut ou si tol=[])
% maxiter : maximum iteration (100 par defaut ou si maxiter=[])
%
% function [u,varargout] = pcg(A,b,tol,maxiter,PC)
% PC : POLYCHAOS sur lequel est defini b et la solution x (A peut etre
% defini sur un autre)
%
% function [u,varargout] = pcg(A,b,tol,maxiter,'precond',precond)
% precond=1 : preconditionneur moyenne (par defaut)
% precond=2 : preconditionneur bloc



if nargin<3 || isempty(tol)
    tol=[];
end
if nargin<4 || isempty(maxiter)
    maxiter =200;
end

if ~israndom(A)
    error('A deterministe : utiliser une autre methode de resolution')
end

[A,b,PC] = pcsystemupdate(A,b,varargin{:});

clock0=clock;

n=size(A,1);

Amean = expect(A);
repzero = find(abs(diag(Amean))/max(abs(diag(Amean)))<1e-12);
repnonzero = setdiff(1:n,repzero);

if ~isempty(repzero)
 fprintf('systeme non defini : resolution sous-systeme\n')
 [utemp,result] = pcg(getcompo(A,repnonzero,repnonzero),getcompo(b,repnonzero,1),tol,maxiter,varargin{:});   
 u = zeros(n,length(PC));
 u(repnonzero,:)=double(utemp);
 u = PCMATRIX(u,[n,1],PC);
    return   
end
masse=getmasse(A);

b=double(cell2mat(PCMATRIX(b)));


precond = getcharin('precond',varargin,1);
switch precond
    case 1
ML = chol(Amean);
    case 2
ML=cell(1,length(PC));
for p=1:length(ML)
ML{p}=sparse(n,n);
for k=1:length(masse)
ML{p}= ML{p} + A.MULTIMATRIX{k}*masse{k}(p,p);
end
ML{p}=chol(ML{p});
end
end

u=zeros(size(b));

if nargout==1
    fprintf('\n');
end


[u(:),flag,relres,iter,resvec]=pcg(@(X) Afun(X,A,n,length(PC)),b(:),tol,maxiter,@(X) Mfun(X,ML,n,length(PC),precond));

if flag~=0
  warning('pcg not converged')  
end
u=PCMATRIX(u,[n,1],PC);

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

