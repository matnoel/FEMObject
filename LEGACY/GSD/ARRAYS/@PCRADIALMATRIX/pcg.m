function [u,result] = pcg(A,b,tol, maxiter,varargin)

if nargin<3 || isempty(tol)
    tol=[];
end
if nargin<4 || isempty(maxiter)
    maxiter = 300;
end

if isa(A,'PCMATRIX')
 [u,result]=pcg(A,PCMATRIX(b),tol,maxiter,varargin{:});  
 return
end

if ~israndom(A)
    error('A deterministe : utiliser une autre methode de resolution')
end



[A,b,PC] = pcsystemupdate(A,b,varargin{:});

clock0=clock;


n=size(A,1);

if isa(b,'double') 
b=b.*one(PC);
end
b=double(expand(cell2mat(b)));


Amean = mean(A);
precond = getcharin('precond',varargin,1);
switch precond
    case 1
ML = chol(Amean,'lower')';
    case 2
ML=cell(1,length(PC));
for p=1:length(ML)
ML{p}=sparse(n,n);
for k=1:length(masse)
    LMk = A.DLmasse{k};
ML{p}= ML{p} + A{k}*LMk(p,p);
end
ML{p}=chol(ML{p},'lower')';
end
end

n=size(A,1);
u=zeros(size(b));



[u(:),flag]=pcg(@(X) Afun(X,A,n,length(PC)),b(:),tol,maxiter,@(X) Mfun(X,ML,n,length(PC),precond));

if flag~=0
  warning('pcg not converged')  
end


u=PCMATRIX(u,[n,1],PC);



if nargout>1
    result.totaltime = etime(clock,clock0);
end



function Y = Afun(X,A,n,P)

X=reshape(full(X),n,P);

temp = X*A.DLmasse;
if iscell(A.V)
temp = mat2cell(temp);
end
Y = multisum(A.V*temp);
Y=full(Y(:));

return

function Y = Mfun(X,ML,n,P,precond)

X = reshape(full(X),n,P);

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

