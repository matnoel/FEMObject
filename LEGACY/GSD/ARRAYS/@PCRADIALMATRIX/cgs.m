function [u,result] = cgs(A,b,tol,maxiter,varargin)
% function [u,result] = cgs(A,b,tol,maxiter,varargin)

if nargin<3
    tol=[];
end
if nargin<4 || isempty(maxiter)
    maxiter = 1000;
end


if isa(A,'PCMATRIX')
    [u,result]=cgs(A,PCMATRIX(b),tol,maxiter,varargin{:});
    return
end

if ~israndom(A)
    error('A deterministe : utiliser une autre methode de resolution')
end

clock0=clock;


Amean = mean(A);
perm=symrcm(Amean);
[ML,MU] = lu(Amean(perm,perm));


[A,b,PC]=pcsystemupdate(A,b,varargin{:});

n=size(A,1);

if isa(b,'double')
    b=b.*one(PC);
end
b=double(cell2mat(PCMATRIX(b)));

u=zeros(size(b));


[u(:),flag]=cgs(@(X) Afun(X,A,n,length(PC)),b(:),tol,maxiter,@(X) Mfun(X,ML,MU,n,length(PC),perm));

if flag~=0
    warning('cgs not converged')
end


u=PCMATRIX(u,[n,1],PC);

% if nargout>1
result.totaltime = etime(clock,clock0);
% end


function Y = Afun(X,A,n,P)


X=reshape(X,n,P);
Y=zeros(n,P);

temp = X*A.DLmasse;
if iscell(A.V)
    temp = mat2cell(temp);
end

Y = multisum(A.V*temp);
Y=Y(:);



return

function Y = Mfun(X,ML,MU,n,P,perm)
Y=reshape(X,n,P);
Y=MU\(ML\Y(perm,:));
Y(perm,:)=Y;
Y = Y(:);
return

